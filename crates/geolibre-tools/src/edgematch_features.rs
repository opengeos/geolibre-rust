//! GeoLibre tool: edgematch line features across a tile/sheet boundary.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Edgematch Links* +
//! *Edgematch Features* (Editing). Completes the conflation suite: GeoLibre
//! ships `integrate`, `rubbersheet_features`, and `detect_feature_changes`, and
//! the bundle provides `transfer_attributes` and `snap_endnodes` — but
//! `snap_endnodes` snaps any nearby endpoints blindly, with no cross-feature
//! one-to-one matching and no notion of a boundary zone. Edgematching (roads or
//! streams digitized per map sheet, meeting at neatlines) needs dangling
//! endpoints paired one-to-one before anything moves.
//!
//! Dangling endpoints (line ends not shared with another feature) within
//! `tolerance` of each other are matched one-to-one by distance (optionally
//! disambiguated by attribute similarity on `match_fields`). Each matched pair is
//! then reconciled: `midpoint` moves both ends to their midpoint, `move_endpoint`
//! snaps the second onto the first. Adjusted lines are written, plus an optional
//! `links` layer of match segments for QA.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct EdgematchFeaturesTool;

impl Tool for EdgematchFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "edgematch_features",
            display_name: "Edgematch Features",
            summary: "Connect line datasets across a tile/sheet boundary (like ArcGIS Generate Edgematch Links + Edgematch Features): match dangling endpoints one-to-one within a tolerance and reconcile them (midpoint or snap), with an optional match-links layer. The cross-feature one-to-one matching the blind snap_endnodes lacks.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line layer (e.g. two adjacent tiles merged) to edgematch.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output adjusted line layer. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Maximum distance between dangling endpoints to match them.",
                    required: true,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'midpoint' (move both ends to their midpoint; default) or 'move_endpoint' (snap the second onto the first).",
                    required: false,
                },
                ToolParamSpec {
                    name: "match_fields",
                    description: "Optional comma-separated fields; a matched pair must agree on these (disambiguates candidates).",
                    required: false,
                },
                ToolParamSpec {
                    name: "links",
                    description: "Optional output line layer of match links (for QA).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let links_path = parse_optional_str(args, "links")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let match_idx: Vec<usize> = prm
            .match_fields
            .iter()
            .filter_map(|f| layer.schema.field_index(f))
            .collect();

        // Collect line endpoints (feature index, end 0=start/1=end, coord).
        let mut ends: Vec<Endpoint> = Vec::new();
        for (fi, feat) in layer.features.iter().enumerate() {
            if let Some(Geometry::LineString(cs)) = feat.geometry.as_ref() {
                if cs.len() >= 2 {
                    ends.push(Endpoint {
                        fi,
                        end: 0,
                        x: cs[0].x,
                        y: cs[0].y,
                    });
                    let last = cs.len() - 1;
                    ends.push(Endpoint {
                        fi,
                        end: 1,
                        x: cs[last].x,
                        y: cs[last].y,
                    });
                }
            }
        }

        // Node valence: an endpoint is dangling if its exact position is not
        // shared by another endpoint.
        let mut node_count: HashMap<(u64, u64), usize> = HashMap::new();
        for e in &ends {
            *node_count
                .entry((e.x.to_bits(), e.y.to_bits()))
                .or_insert(0) += 1;
        }
        let dangling: Vec<usize> = (0..ends.len())
            .filter(|&i| node_count[&(ends[i].x.to_bits(), ends[i].y.to_bits())] == 1)
            .collect();

        ctx.progress
            .info(&format!("{} dangling endpoint(s) to match", dangling.len()));

        // Candidate pairs: different features, within tolerance, distance > 0,
        // agreeing on match_fields. Greedy one-to-one by ascending distance.
        let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
        for a in 0..dangling.len() {
            for b in (a + 1)..dangling.len() {
                let ea = &ends[dangling[a]];
                let eb = &ends[dangling[b]];
                if ea.fi == eb.fi {
                    continue;
                }
                let d = (ea.x - eb.x).hypot(ea.y - eb.y);
                if d <= 0.0 || d > prm.tolerance {
                    continue;
                }
                if !attrs_agree(&layer, ea.fi, eb.fi, &match_idx) {
                    continue;
                }
                pairs.push((d, dangling[a], dangling[b]));
            }
        }
        pairs.sort_by(|x, y| x.0.total_cmp(&y.0));

        let mut used = vec![false; ends.len()];
        let mut matches: Vec<(usize, usize)> = Vec::new();
        for (_, a, b) in pairs {
            if used[a] || used[b] {
                continue;
            }
            used[a] = true;
            used[b] = true;
            matches.push((a, b));
        }

        // Apply adjustments to the geometry.
        let mut link_layer = Layer::new("links").with_geom_type(GeometryType::LineString);
        if let Some(e) = layer.crs_epsg() {
            link_layer = link_layer.with_crs_epsg(e);
        }
        link_layer.add_field(FieldDef::new("dist", FieldType::Float));

        for &(a, b) in &matches {
            let (ax, ay) = (ends[a].x, ends[a].y);
            let (bx, by) = (ends[b].x, ends[b].y);
            let (nx, ny) = match prm.method {
                Method::Midpoint => ((ax + bx) / 2.0, (ay + by) / 2.0),
                Method::MoveEndpoint => (ax, ay),
            };
            set_endpoint(&mut layer.features[ends[a].fi], ends[a].end, nx, ny);
            set_endpoint(&mut layer.features[ends[b].fi], ends[b].end, nx, ny);
            if links_path.is_some() {
                link_layer.push(Feature {
                    fid: 0,
                    geometry: Some(Geometry::line_string(vec![
                        Coord::xy(ax, ay),
                        Coord::xy(bx, by),
                    ])),
                    attributes: vec![FieldValue::Float((ax - bx).hypot(ay - by))],
                });
            }
        }

        let matched = matches.len();
        let out_path = write_or_store_layer(layer, output)?;
        let links_out = match links_path {
            Some(p) => Some(write_or_store_layer(link_layer, Some(p))?),
            None => None,
        };

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("dangling_endpoints".to_string(), json!(dangling.len()));
        outputs.insert("matches".to_string(), json!(matched));
        if let Some(l) = links_out {
            outputs.insert("links".to_string(), json!(l));
        }
        Ok(ToolRunResult { outputs })
    }
}

fn set_endpoint(feat: &mut Feature, end: u8, x: f64, y: f64) {
    if let Some(Geometry::LineString(cs)) = feat.geometry.as_mut() {
        let idx = if end == 0 { 0 } else { cs.len() - 1 };
        cs[idx] = Coord::xy(x, y);
    }
}

fn attrs_agree(layer: &Layer, fa: usize, fb: usize, idx: &[usize]) -> bool {
    for &i in idx {
        let va = layer.features[fa].attributes.get(i);
        let vb = layer.features[fb].attributes.get(i);
        if value_string_opt(va) != value_string_opt(vb) {
            return false;
        }
    }
    true
}

fn value_string_opt(fv: Option<&FieldValue>) -> String {
    match fv {
        Some(v) if v.as_i64().is_some() => v.as_i64().unwrap().to_string(),
        Some(v) if v.as_f64().is_some() => format!("{}", v.as_f64().unwrap()),
        Some(v) => v.as_str().unwrap_or("").to_string(),
        None => String::new(),
    }
}

struct Endpoint {
    fi: usize,
    end: u8,
    x: f64,
    y: f64,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Midpoint,
    MoveEndpoint,
}

struct Params {
    tolerance: f64,
    method: Method,
    match_fields: Vec<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let tolerance = match args.get("tolerance") {
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'tolerance' must be a number".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "required parameter 'tolerance' is missing".into(),
            ))
        }
    };
    if !(tolerance > 0.0) {
        return Err(ToolError::Validation("'tolerance' must be positive".into()));
    }
    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("midpoint") => Method::Midpoint,
        Some("move_endpoint") => Method::MoveEndpoint,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'midpoint' or 'move_endpoint', got '{o}'"
            )))
        }
    };
    let match_fields = match args.get("match_fields").and_then(Value::as_str) {
        None => Vec::new(),
        Some(s) => s
            .split(',')
            .map(str::trim)
            .filter(|x| !x.is_empty())
            .map(String::from)
            .collect(),
    };
    Ok(Params {
        tolerance,
        method,
        match_fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line_layer(lines: &[(&str, Vec<(f64, f64)>)]) -> String {
        let mut l = Layer::new("ln")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("road", FieldType::Text));
        for (road, pts) in lines {
            let coords: Vec<Coord> = pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(
                Some(Geometry::line_string(coords)),
                &[("road", (*road).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EdgematchFeaturesTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// Two lines whose ends nearly meet across a boundary get matched and joined
    /// at the midpoint.
    #[test]
    fn matches_and_joins_endpoints() {
        // Line A ends at (10, 0.4); line B starts at (10, -0.4). Gap 0.8.
        let input = line_layer(&[
            ("main", vec![(0.0, 0.0), (10.0, 0.4)]),
            ("main", vec![(10.0, -0.4), (20.0, 0.0)]),
        ]);
        let (out, l) = run(json!({ "input": input, "tolerance": 2.0, "method": "midpoint" }));
        assert_eq!(
            out.outputs["matches"],
            json!(1),
            "the two dangling ends should match"
        );
        // Both endpoints now at the midpoint (10, 0).
        let a_end = last_pt(&l.features[0]);
        let b_start = first_pt(&l.features[1]);
        assert!((a_end.0 - 10.0).abs() < 1e-9 && (a_end.1 - 0.0).abs() < 1e-9);
        assert_eq!(a_end, b_start, "the ends must now coincide");
    }

    /// Ends farther apart than the tolerance are not matched.
    #[test]
    fn respects_tolerance() {
        let input = line_layer(&[
            ("a", vec![(0.0, 0.0), (10.0, 5.0)]),
            ("b", vec![(10.0, -5.0), (20.0, 0.0)]),
        ]);
        let (out, _l) = run(json!({ "input": input, "tolerance": 2.0 }));
        assert_eq!(out.outputs["matches"], json!(0), "gap 10 > tolerance 2");
    }

    /// match_fields prevents joining lines of different classes.
    #[test]
    fn match_fields_disambiguate() {
        let input = line_layer(&[
            ("hwy", vec![(0.0, 0.0), (10.0, 0.3)]),
            ("local", vec![(10.0, -0.3), (20.0, 0.0)]),
        ]);
        let (out, _l) = run(json!({
            "input": input, "tolerance": 2.0, "match_fields": "road",
        }));
        assert_eq!(
            out.outputs["matches"],
            json!(0),
            "different road classes must not match"
        );
    }

    fn first_pt(f: &Feature) -> (f64, f64) {
        if let Some(Geometry::LineString(cs)) = &f.geometry {
            (cs[0].x, cs[0].y)
        } else {
            unreachable!()
        }
    }
    fn last_pt(f: &Feature) -> (f64, f64) {
        if let Some(Geometry::LineString(cs)) = &f.geometry {
            (cs[cs.len() - 1].x, cs[cs.len() - 1].y)
        } else {
            unreachable!()
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            EdgematchFeaturesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no tolerance
        assert!(bad(json!({ "input": "a.geojson", "tolerance": 0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "tolerance": 5, "method": "weld" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "tolerance": 5 })).is_ok());
    }
}
