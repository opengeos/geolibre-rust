//! GeoLibre tool: origin-destination "desire lines".
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Origin-Destination Links*
//! (Analysis). The bundled OD tools (`network_od_cost_matrix`,
//! `multimodal_od_cost_matrix`) output cost *tables*, not geometry; nothing
//! draws the straight spider/desire lines used for flow mapping and catchment
//! visualization. This connects origin points to destination points and emits
//! one `LineString` per pair, each carrying the origin id, destination id, and
//! length — ready for `render_vector_png` or `vector_to_pmtiles`.
//!
//! Pairing is chosen by three combinable rules:
//!
//! * `id_field` — link only origin/destination pairs whose shared field matches
//!   (e.g. link each trip to the store it was assigned to);
//! * `search_distance` — keep only destinations within a radius;
//! * `num_nearest` — keep the k nearest of the surviving destinations.
//!
//! With none set, each origin links to its single nearest destination.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GenerateOdLinksTool;

impl Tool for GenerateOdLinksTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_od_links",
            display_name: "Generate Origin-Destination Links",
            summary: "Draw straight desire lines from origin points to destination points — matched by a shared id, within a distance, and/or to the k nearest — each carrying origin/destination ids and length, like ArcGIS Generate Origin-Destination Links.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "origins",
                    description: "Origin point layer (Point or MultiPoint).",
                    required: true,
                },
                ToolParamSpec {
                    name: "destinations",
                    description: "Destination point layer (Point or MultiPoint).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_nearest",
                    description: "Keep only the k nearest destinations per origin (after any id/distance filter).",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Keep only destinations within this distance of the origin (map units).",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Field present in BOTH layers; link only pairs whose values match. Also written as 'link_id'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "origin_id_field",
                    description: "Field in origins to write as 'origin_id' (default: the feature index).",
                    required: false,
                },
                ToolParamSpec {
                    name: "dest_id_field",
                    description: "Field in destinations to write as 'dest_id' (default: the feature index).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["origins", "destinations"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let origins_path = args.get("origins").and_then(Value::as_str).unwrap();
        let dests_path = args.get("destinations").and_then(Value::as_str).unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let origins = load_input_layer(origins_path)?;
        let dests = load_input_layer(dests_path)?;

        let o_pts = collect_points(
            &origins,
            prm.origin_id_field.as_deref(),
            prm.id_field.as_deref(),
        )?;
        let d_pts = collect_points(
            &dests,
            prm.dest_id_field.as_deref(),
            prm.id_field.as_deref(),
        )?;
        if o_pts.is_empty() || d_pts.is_empty() {
            return Err(ToolError::Execution(
                "both origins and destinations must contain point features".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "linking {} origins to {} destinations",
            o_pts.len(),
            d_pts.len()
        ));

        // Default: nearest destination only when NO pairing rule is given.
        let default_nearest =
            prm.num_nearest.is_none() && prm.search_distance.is_none() && prm.id_field.is_none();

        let mut out = Layer::new("od_links").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = origins.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("origin_id", FieldType::Text));
        out.add_field(FieldDef::new("dest_id", FieldType::Text));
        out.add_field(FieldDef::new("length", FieldType::Float));
        if prm.id_field.is_some() {
            out.add_field(FieldDef::new("link_id", FieldType::Text));
        }

        let radius2 = prm.search_distance.map(|d| d * d);
        let mut link_count = 0usize;
        let mut total_len = 0.0;
        for o in &o_pts {
            // Candidate destinations: match id, then distance filter.
            let mut cands: Vec<(f64, &Pt)> = Vec::new();
            for d in &d_pts {
                if let (Some(oid), Some(did)) = (&o.link, &d.link) {
                    if oid != did {
                        continue;
                    }
                } else if prm.id_field.is_some() {
                    // id_field requested but one side lacks a value -> no match.
                    continue;
                }
                let dist2 = (o.x - d.x).powi(2) + (o.y - d.y).powi(2);
                if let Some(r2) = radius2 {
                    if dist2 > r2 {
                        continue;
                    }
                }
                cands.push((dist2, d));
            }

            // Keep k nearest (or 1 by default) if requested.
            let k = if let Some(k) = prm.num_nearest {
                Some(k)
            } else if default_nearest {
                Some(1)
            } else {
                None
            };
            if let Some(k) = k {
                cands.sort_by(|a, b| a.0.total_cmp(&b.0));
                cands.truncate(k);
            }

            for (dist2, d) in cands {
                let length = dist2.sqrt();
                let geom = Geometry::LineString(vec![Coord::xy(o.x, o.y), Coord::xy(d.x, d.y)]);
                let mut attrs: Vec<(&str, FieldValue)> = vec![
                    ("origin_id", FieldValue::Text(o.id.clone())),
                    ("dest_id", FieldValue::Text(d.id.clone())),
                    ("length", FieldValue::Float(length)),
                ];
                let link_val;
                if prm.id_field.is_some() {
                    link_val = o.link.clone().unwrap_or_default();
                    attrs.push(("link_id", FieldValue::Text(link_val)));
                }
                out.add_feature(Some(geom), &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed writing link: {e}")))?;
                link_count += 1;
                total_len += length;
            }
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("origin_count".to_string(), json!(o_pts.len()));
        outputs.insert("destination_count".to_string(), json!(d_pts.len()));
        outputs.insert("link_count".to_string(), json!(link_count));
        outputs.insert(
            "mean_length".to_string(),
            json!(if link_count > 0 {
                total_len / link_count as f64
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

struct Pt {
    x: f64,
    y: f64,
    id: String,           // label for output
    link: Option<String>, // matching value (id_field), if any
}

/// Collect points with an id label (`id_field` for output) and an optional
/// match value (`link_field`, shared across both layers).
fn collect_points(
    layer: &Layer,
    id_field: Option<&str>,
    link_field: Option<&str>,
) -> Result<Vec<Pt>, ToolError> {
    let id_idx = match id_field {
        Some(f) => Some(
            layer
                .schema
                .field_index(f)
                .ok_or_else(|| ToolError::Validation(format!("id field '{f}' not found")))?,
        ),
        None => None,
    };
    let link_idx = match link_field {
        Some(f) => layer.schema.field_index(f), // absent -> None (no match on this side)
        None => None,
    };

    let mut pts = Vec::new();
    for (fidx, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        let coords: Vec<(f64, f64)> = match geom {
            Geometry::Point(c) => vec![(c.x, c.y)],
            Geometry::MultiPoint(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
            _ => continue,
        };
        let id = match id_idx {
            Some(i) => field_key(&feature.attributes[i]),
            None => fidx.to_string(),
        };
        let link = link_idx.map(|i| field_key(&feature.attributes[i]));
        for (x, y) in coords {
            pts.push(Pt {
                x,
                y,
                id: id.clone(),
                link: link.clone(),
            });
        }
    }
    Ok(pts)
}

fn field_key(fv: &FieldValue) -> String {
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
    num_nearest: Option<usize>,
    search_distance: Option<f64>,
    id_field: Option<String>,
    origin_id_field: Option<String>,
    dest_id_field: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let num_nearest = match args.get("num_nearest") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => Some(n.as_u64().unwrap_or(0).max(1) as usize),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<usize>()
                .map_err(|_| ToolError::Validation("'num_nearest' must be an integer".into()))?
                .max(1),
        ),
        Some(_) => {
            return Err(ToolError::Validation(
                "'num_nearest' must be a number".to_string(),
            ))
        }
    };
    let search_distance = match args.get("search_distance") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_f64(),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<f64>()
                .map_err(|_| ToolError::Validation("'search_distance' must be a number".into()))?,
        ),
        Some(_) => {
            return Err(ToolError::Validation(
                "'search_distance' must be a number".to_string(),
            ))
        }
    };
    if let Some(d) = search_distance {
        if d.is_nan() || d <= 0.0 {
            return Err(ToolError::Validation(
                "'search_distance' must be positive".to_string(),
            ));
        }
    }
    Ok(Params {
        num_nearest,
        search_distance,
        id_field: parse_optional_str(args, "id_field")?.map(str::to_string),
        origin_id_field: parse_optional_str(args, "origin_id_field")?.map(str::to_string),
        dest_id_field: parse_optional_str(args, "dest_id_field")?.map(str::to_string),
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

    fn point_layer(pts: &[(f64, f64, &str)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("zone", FieldType::Text));
        for (x, y, z) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("zone", (*z).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GenerateOdLinksTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Default (no rule): each origin links to its single nearest destination.
    #[test]
    fn default_links_to_nearest() {
        let o = point_layer(&[(0.0, 0.0, "a")]);
        let d = point_layer(&[(3.0, 4.0, "x"), (100.0, 0.0, "y")]);
        let (out, layer) = run(json!({ "origins": o, "destinations": d }));
        assert_eq!(out.outputs["link_count"], json!(1));
        let lidx = layer.schema.field_index("length").unwrap();
        let f = layer.iter().next().unwrap();
        assert!(
            (f.attributes[lidx].as_f64().unwrap() - 5.0).abs() < 1e-6,
            "3-4-5"
        );
    }

    /// num_nearest = k links each origin to its k closest destinations.
    #[test]
    fn k_nearest_links() {
        let o = point_layer(&[(0.0, 0.0, "a")]);
        let d = point_layer(&[(1.0, 0.0, "x"), (2.0, 0.0, "y"), (50.0, 0.0, "z")]);
        let (out, _l) = run(json!({ "origins": o, "destinations": d, "num_nearest": 2 }));
        assert_eq!(out.outputs["link_count"], json!(2));
    }

    /// search_distance keeps only destinations within the radius.
    #[test]
    fn radius_filters_links() {
        let o = point_layer(&[(0.0, 0.0, "a")]);
        let d = point_layer(&[(3.0, 0.0, "x"), (8.0, 0.0, "y"), (50.0, 0.0, "z")]);
        let (out, _l) = run(json!({ "origins": o, "destinations": d, "search_distance": 10.0 }));
        // 3 and 8 within 10; 50 excluded. No num_nearest -> all within radius.
        assert_eq!(out.outputs["link_count"], json!(2));
    }

    /// id_field links only pairs whose shared field matches.
    #[test]
    fn id_field_matches_pairs() {
        let o = point_layer(&[(0.0, 0.0, "north"), (0.0, 100.0, "south")]);
        let d = point_layer(&[
            (10.0, 0.0, "north"),
            (10.0, 100.0, "south"),
            (5.0, 5.0, "north"),
        ]);
        let (out, layer) = run(json!({ "origins": o, "destinations": d, "id_field": "zone" }));
        // north origin -> 2 north dests; south origin -> 1 south dest = 3 links.
        assert_eq!(out.outputs["link_count"], json!(3));
        let idx = layer.schema.field_index("link_id").unwrap();
        for f in layer.iter() {
            let v = f.attributes[idx].as_str().unwrap();
            assert!(v == "north" || v == "south");
        }
    }

    /// id_field combined with num_nearest keeps the nearest matching dest.
    #[test]
    fn id_field_with_num_nearest() {
        let o = point_layer(&[(0.0, 0.0, "north")]);
        let d = point_layer(&[
            (10.0, 0.0, "north"),
            (5.0, 0.0, "north"),
            (1.0, 0.0, "south"),
        ]);
        let (out, layer) =
            run(json!({ "origins": o, "destinations": d, "id_field": "zone", "num_nearest": 1 }));
        assert_eq!(out.outputs["link_count"], json!(1));
        let lidx = layer.schema.field_index("length").unwrap();
        let f = layer.iter().next().unwrap();
        // nearest NORTH dest is at x=5 (south x=1 excluded by id match).
        assert!((f.attributes[lidx].as_f64().unwrap() - 5.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_missing_destinations() {
        let o = point_layer(&[(0.0, 0.0, "a")]);
        let args: ToolArgs = serde_json::from_value(json!({ "origins": o })).unwrap();
        assert!(GenerateOdLinksTool.validate(&args).is_err());
    }
}
