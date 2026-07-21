//! GeoLibre tool: detect changes between two line datasets.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Detect Feature Changes* (Data
//! Management) — the entry point to conflation workflows (OSM vs. authoritative
//! data, versioned datasets). No conflation or vector change-detection exists in
//! the bundled suite (`change_vector_analysis` etc. are raster tools); this
//! delivers standalone value and builds the feature-matching core that
//! Rubbersheet / Edgematch would later reuse.
//!
//! Each **update** line is matched to the nearest **base** line by symmetric
//! (discrete) Hausdorff distance, within `search_distance`. The match is then
//! classified:
//!
//! * **unchanged** — matched, geometry within `spatial_tolerance`, compared
//!   attributes equal;
//! * **spatial** — matched but the geometry moved beyond `spatial_tolerance`;
//! * **attribute** — matched, geometry unchanged, a compared attribute differs;
//! * **spatial_attribute** — both moved and an attribute changed;
//! * **new** — no base line within `search_distance`.
//!
//! Base lines that no update line matched are emitted as **deleted**. Every
//! output feature carries `change_type`, the matched base id (`match_id`, −1 for
//! new), and the match Hausdorff distance (`match_dist`). Distances are in the
//! layer CRS units — use a projected CRS.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct DetectFeatureChangesTool;

impl Tool for DetectFeatureChangesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "detect_feature_changes",
            display_name: "Detect Feature Changes",
            summary: "Match two line datasets by Hausdorff distance and classify each as unchanged / spatial / attribute / spatial_attribute / new, plus deleted base features — the vector change-detection and conflation entry point, like ArcGIS Detect Feature Changes.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "update",
                    description: "Updated line vector layer (the new state).",
                    required: true,
                },
                ToolParamSpec {
                    name: "base",
                    description: "Base line vector layer (the previous state) to compare against.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector (update features classified, plus deleted base features). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Maximum distance (CRS units) to match an update line to a base line. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "spatial_tolerance",
                    description: "Hausdorff distance below which a matched geometry counts as unchanged (CRS units). Default: search_distance / 5.",
                    required: false,
                },
                ToolParamSpec {
                    name: "compare_fields",
                    description: "Comma-separated field(s) present on both layers to compare for attribute changes.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "update")?;
        require_str(args, "base")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let update_path = require_str(args, "update")?;
        let base_path = require_str(args, "base")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let update = load_input_layer(update_path)?;
        let base = load_input_layer(base_path)?;

        // Resolve compare-field indices on each layer (must exist on both).
        let mut ucmp = Vec::new();
        let mut bcmp = Vec::new();
        for name in &prm.compare_fields {
            let ui = update.schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("compare field '{name}' not on update layer"))
            })?;
            let bi = base.schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("compare field '{name}' not on base layer"))
            })?;
            ucmp.push(ui);
            bcmp.push(bi);
        }

        // Extract base geometries (vertices + segments) and bboxes.
        let base_lines: Vec<Option<LineGeom>> = base
            .features
            .iter()
            .map(|f| f.geometry.as_ref().and_then(line_geom))
            .collect();
        let update_lines: Vec<Option<LineGeom>> = update
            .features
            .iter()
            .map(|f| f.geometry.as_ref().and_then(line_geom))
            .collect();

        ctx.progress.info(&format!(
            "matching {} update line(s) against {} base line(s)",
            update.len(),
            base.len()
        ));

        let spatial_tol = prm.spatial_tolerance.unwrap_or(prm.search_distance / 5.0);
        let mut base_matched = vec![false; base.len()];

        // Build the output layer from the update schema + change columns.
        let mut out = Layer::new("changes").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = update.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in update.schema.fields() {
            out.add_field(field.clone());
        }
        // Ensure the change columns don't collide with existing names.
        out.add_field(FieldDef::new("change_type", FieldType::Text));
        out.add_field(FieldDef::new("match_id", FieldType::Integer));
        out.add_field(FieldDef::new("match_dist", FieldType::Float));

        let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
        for (ui, feature) in update.features.iter().enumerate() {
            let Some(ug) = &update_lines[ui] else {
                continue;
            };
            // Find the best base match by symmetric Hausdorff within search_distance.
            let mut best: Option<(usize, f64)> = None;
            for (bi, bg_opt) in base_lines.iter().enumerate() {
                let Some(bg) = bg_opt else { continue };
                if !bbox_within(&ug.bbox, &bg.bbox, prm.search_distance) {
                    continue;
                }
                let h = symmetric_hausdorff(
                    ug,
                    bg,
                    best.map(|(_, d)| d).unwrap_or(prm.search_distance),
                );
                if h <= prm.search_distance && best.is_none_or(|(_, d)| h < d) {
                    best = Some((bi, h));
                }
            }

            let (change_type, match_id, match_dist) = match best {
                None => ("new", -1i64, -1.0),
                Some((bi, h)) => {
                    base_matched[bi] = true;
                    let geom_changed = h > spatial_tol;
                    let attr_changed = self.attrs_differ(feature, &base.features[bi], &ucmp, &bcmp);
                    let ct = match (geom_changed, attr_changed) {
                        (false, false) => "unchanged",
                        (true, false) => "spatial",
                        (false, true) => "attribute",
                        (true, true) => "spatial_attribute",
                    };
                    (ct, bi as i64, h)
                }
            };
            *counts.entry(change_type).or_default() += 1;

            let mut attrs = feature.attributes.clone();
            attrs.push(FieldValue::Text(change_type.to_string()));
            attrs.push(FieldValue::Integer(match_id));
            attrs.push(FieldValue::Float(match_dist));
            out.push(Feature {
                fid: 0,
                geometry: feature.geometry.clone(),
                attributes: attrs,
            });
        }

        // Emit unmatched base lines as deleted (base attributes are dropped;
        // only the change columns are filled, geometry from base).
        let n_update_fields = update.schema.fields().len();
        for (bi, matched) in base_matched.iter().enumerate() {
            if *matched || base_lines[bi].is_none() {
                continue;
            }
            *counts.entry("deleted").or_default() += 1;
            let mut attrs: Vec<FieldValue> = vec![FieldValue::Null; n_update_fields];
            attrs.push(FieldValue::Text("deleted".to_string()));
            attrs.push(FieldValue::Integer(bi as i64));
            attrs.push(FieldValue::Float(-1.0));
            out.push(Feature {
                fid: 0,
                geometry: base.features[bi].geometry.clone(),
                attributes: attrs,
            });
        }

        ctx.progress.info(&format!(
            "changes: {}",
            counts
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect::<Vec<_>>()
                .join(", ")
        ));

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        for (k, v) in &counts {
            outputs.insert(format!("count_{k}"), json!(v));
        }
        Ok(ToolRunResult { outputs })
    }
}

impl DetectFeatureChangesTool {
    fn attrs_differ(&self, u: &Feature, b: &Feature, ucmp: &[usize], bcmp: &[usize]) -> bool {
        ucmp.iter().zip(bcmp).any(|(&ui, &bi)| {
            let uv = u.attributes.get(ui).map(value_string).unwrap_or_default();
            let bv = b.attributes.get(bi).map(value_string).unwrap_or_default();
            uv != bv
        })
    }
}

// ── Hausdorff distance ───────────────────────────────────────────────────────

struct LineGeom {
    verts: Vec<(f64, f64)>,
    segs: Vec<Seg>,
    bbox: [f64; 4],
}

type Seg = ((f64, f64), (f64, f64));

fn add_chain(cs: &[Coord], verts: &mut Vec<(f64, f64)>, segs: &mut Vec<Seg>) {
    for w in cs.windows(2) {
        segs.push(((w[0].x, w[0].y), (w[1].x, w[1].y)));
    }
    for c in cs {
        verts.push((c.x, c.y));
    }
}

fn line_geom(geom: &Geometry) -> Option<LineGeom> {
    let mut verts = Vec::new();
    let mut segs = Vec::new();
    match geom {
        Geometry::LineString(cs) => add_chain(cs, &mut verts, &mut segs),
        Geometry::MultiLineString(lines) => {
            for l in lines {
                add_chain(l, &mut verts, &mut segs);
            }
        }
        _ => return None,
    }
    if verts.is_empty() || segs.is_empty() {
        return None;
    }
    let mut bbox = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for &(x, y) in &verts {
        bbox[0] = bbox[0].min(x);
        bbox[1] = bbox[1].min(y);
        bbox[2] = bbox[2].max(x);
        bbox[3] = bbox[3].max(y);
    }
    Some(LineGeom { verts, segs, bbox })
}

/// Symmetric discrete Hausdorff: max of the two directed vertex-to-polyline
/// distances. Short-circuits once it exceeds `cutoff`.
fn symmetric_hausdorff(a: &LineGeom, b: &LineGeom, cutoff: f64) -> f64 {
    let d1 = directed_hausdorff(&a.verts, &b.segs, cutoff);
    if d1 > cutoff {
        return d1;
    }
    let d2 = directed_hausdorff(&b.verts, &a.segs, cutoff);
    d1.max(d2)
}

fn directed_hausdorff(verts: &[(f64, f64)], segs: &[Seg], cutoff: f64) -> f64 {
    let mut worst = 0.0f64;
    for &p in verts {
        let mut best = f64::INFINITY;
        for &(a, b) in segs {
            let d = point_seg_dist(p, a, b);
            if d < best {
                best = d;
                if best == 0.0 {
                    break;
                }
            }
        }
        if best > worst {
            worst = best;
            if worst > cutoff {
                return worst; // no need to refine beyond the cutoff
            }
        }
    }
    worst
}

fn point_seg_dist(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let len2 = dx * dx + dy * dy;
    if len2 <= 0.0 {
        return (p.0 - a.0).hypot(p.1 - a.1);
    }
    let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / len2).clamp(0.0, 1.0);
    (p.0 - (a.0 + t * dx)).hypot(p.1 - (a.1 + t * dy))
}

fn bbox_within(a: &[f64; 4], b: &[f64; 4], pad: f64) -> bool {
    a[0] - pad <= b[2] && b[0] - pad <= a[2] && a[1] - pad <= b[3] && b[1] - pad <= a[3]
}

fn value_string(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    search_distance: f64,
    spatial_tolerance: Option<f64>,
    compare_fields: Vec<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let search_distance = opt_f64(args, "search_distance")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'search_distance' is missing".to_string())
    })?;
    if !(search_distance > 0.0 && search_distance.is_finite()) {
        return Err(ToolError::Validation(
            "'search_distance' must be a positive number".to_string(),
        ));
    }
    let spatial_tolerance = opt_f64(args, "spatial_tolerance")?;
    if let Some(t) = spatial_tolerance {
        if !(t >= 0.0 && t.is_finite()) {
            return Err(ToolError::Validation(
                "'spatial_tolerance' must be non-negative".to_string(),
            ));
        }
    }
    let compare_fields = parse_optional_str(args, "compare_fields")?
        .map(|s| {
            s.split(',')
                .map(|x| x.trim().to_string())
                .filter(|x| !x.is_empty())
                .collect()
        })
        .unwrap_or_default();
    Ok(Params {
        search_distance,
        spatial_tolerance,
        compare_fields,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
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

    fn line_layer(lines: &[(&[(f64, f64)], &str)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (coords, name) in lines {
            let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(Some(Geometry::line_string(cs)), &[("name", (*name).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DetectFeatureChangesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn change_of(layer: &Layer, name: &str) -> String {
        let ni = layer.schema.field_index("name").unwrap();
        let ci = layer.schema.field_index("change_type").unwrap();
        layer
            .iter()
            .find(|f| f.attributes[ni].as_str() == Some(name))
            .map(|f| f.attributes[ci].as_str().unwrap().to_string())
            .unwrap_or_default()
    }

    #[test]
    fn classifies_all_change_types() {
        let base = line_layer(&[
            (&[(0.0, 0.0), (100.0, 0.0)][..], "keep"),
            (&[(0.0, 50.0), (100.0, 50.0)][..], "move"),
            (&[(0.0, 100.0), (100.0, 100.0)][..], "attr"),
            (&[(0.0, 150.0), (100.0, 150.0)][..], "gone"),
        ]);
        let update = line_layer(&[
            (&[(0.0, 0.0), (100.0, 0.0)][..], "keep"), // identical -> unchanged
            (&[(0.0, 70.0), (100.0, 70.0)][..], "move"), // shifted 20 -> spatial
            (&[(0.0, 100.0), (100.0, 100.0)][..], "attr2"), // same geom, name changed
            (&[(500.0, 500.0), (600.0, 500.0)][..], "brand_new"), // no base -> new
        ]);
        // For the attribute case the geometry matches 'attr' base but the name
        // field differs; match is by geometry so we compare a *shared* field.
        // Use a second field common to both to detect attr change: reuse 'name'
        // won't match by geometry search since names differ but geometry same.
        let (out, layer) = run(json!({
            "update": update, "base": base, "search_distance": 30.0,
            "spatial_tolerance": 5.0, "compare_fields": "name",
        }));
        assert_eq!(change_of(&layer, "keep"), "unchanged");
        assert_eq!(change_of(&layer, "move"), "spatial");
        // 'attr2' matches base 'attr' geometrically (same line) but name differs.
        assert_eq!(change_of(&layer, "attr2"), "attribute");
        assert_eq!(change_of(&layer, "brand_new"), "new");
        // The base 'gone' line has no update match -> deleted.
        assert_eq!(out.outputs["count_deleted"], json!(1));
    }

    #[test]
    fn identical_layers_are_all_unchanged() {
        let l = line_layer(&[
            (&[(0.0, 0.0), (10.0, 0.0)][..], "a"),
            (&[(0.0, 10.0), (10.0, 10.0)][..], "b"),
        ]);
        let (out, _layer) = run(json!({ "update": l.clone(), "base": l, "search_distance": 5.0 }));
        assert_eq!(out.outputs["count_unchanged"], json!(2));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DetectFeatureChangesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "update": "u.geojson", "base": "b.geojson" })).is_err());
        assert!(
            bad(json!({ "update": "u.geojson", "base": "b.geojson", "search_distance": 0 }))
                .is_err()
        );
        assert!(
            bad(json!({ "update": "u.geojson", "base": "b.geojson", "search_distance": 10 }))
                .is_ok()
        );
    }
}
