//! GeoLibre tool: `adjust_3d_z` — apply a Z-unit conversion, vertical offset,
//! and/or multiplicative factor to the Z ordinates of already-3D features.
//!
//! ArcGIS counterpart: **Adjust 3D Z**
//! <https://pro.arcgis.com/en/pro-app/latest/tool-reference/data-management/adjust-3d-z.htm>
//!
//! The transform applied to every Z ordinate is:
//!
//! ```text
//! z' = z * unit_conversion * factor + offset
//! ```
//!
//! where `unit_conversion` is derived from the `from_unit`/`to_unit` presets
//! (e.g. feet → meters = 0.3048), `factor` is an explicit multiplier (vertical
//! exaggeration or a datum scale), and `offset` is an additive vertical shift
//! expressed in the target units. Vertices without a Z value (2D features, or 2D
//! vertices of a mixed geometry) pass through unchanged.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Applies a Z-unit conversion, multiplicative factor, and vertical offset to
/// the Z ordinates of 3D vector features.
pub struct Adjust3dZTool;

impl Tool for Adjust3dZTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "adjust_3d_z",
            display_name: "Adjust 3D Z",
            summary: "Apply a Z-unit conversion, vertical offset, and/or multiplicative factor to the Z ordinates of 3D features.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (features should carry Z ordinates).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "factor",
                    description: "Multiplicative factor applied to Z (default 1.0). Useful for vertical exaggeration or a datum scale.",
                    required: false,
                },
                ToolParamSpec {
                    name: "offset",
                    description: "Additive vertical offset applied after the factor, in target units (default 0.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "from_unit",
                    description: "Source Z unit preset (meters, feet, us_feet, centimeters, millimeters, kilometers, miles, yards, inches). Requires to_unit.",
                    required: false,
                },
                ToolParamSpec {
                    name: "to_unit",
                    description: "Target Z unit preset. Requires from_unit. Z is scaled by (meters-per-from_unit / meters-per-to_unit).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        // Parse all optional params up front so bad values fail fast.
        Params::parse(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_str(args, "output")?;
        let params = Params::parse(args)?;

        let combined_factor = params.unit_conversion * params.factor;
        let offset = params.offset;

        ctx.progress.info(&format!(
            "adjusting Z: z' = z * {combined_factor} + {offset}"
        ));

        let mut layer = load_input_layer(input)?;

        let mut features_with_z: u64 = 0;
        let mut vertices_adjusted: u64 = 0;
        let mut z_min = f64::INFINITY;
        let mut z_max = f64::NEG_INFINITY;

        let n = layer.features.len().max(1);
        for (i, feature) in layer.features.iter_mut().enumerate() {
            if let Some(geom) = feature.geometry.as_mut() {
                let mut feature_has_z = false;
                for_each_coord_mut(geom, &mut |c: &mut Coord| {
                    if let Some(z) = c.z {
                        let nz = z * combined_factor + offset;
                        c.z = Some(nz);
                        feature_has_z = true;
                        vertices_adjusted += 1;
                        if nz < z_min {
                            z_min = nz;
                        }
                        if nz > z_max {
                            z_max = nz;
                        }
                    }
                });
                if feature_has_z {
                    features_with_z += 1;
                }
            }
            ctx.progress.progress((i as f64 + 1.0) / n as f64);
        }

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("features_with_z".to_string(), json!(features_with_z));
        outputs.insert("vertices_adjusted".to_string(), json!(vertices_adjusted));
        outputs.insert("combined_factor".to_string(), json!(combined_factor));
        outputs.insert("offset".to_string(), json!(offset));
        if vertices_adjusted > 0 {
            outputs.insert("z_min".to_string(), json!(z_min));
            outputs.insert("z_max".to_string(), json!(z_max));
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Parsed / validated parameters for a run.
struct Params {
    factor: f64,
    offset: f64,
    unit_conversion: f64,
}

impl Params {
    fn parse(args: &ToolArgs) -> Result<Self, ToolError> {
        let factor = parse_optional_f64(args, "factor")?.unwrap_or(1.0);
        let offset = parse_optional_f64(args, "offset")?.unwrap_or(0.0);

        let from_unit = parse_optional_str(args, "from_unit")?;
        let to_unit = parse_optional_str(args, "to_unit")?;

        let unit_conversion = match (from_unit, to_unit) {
            (None, None) => 1.0,
            (Some(f), Some(t)) => {
                let mf = meters_per_unit(f)?;
                let mt = meters_per_unit(t)?;
                mf / mt
            }
            (Some(_), None) | (None, Some(_)) => {
                return Err(ToolError::Validation(
                    "'from_unit' and 'to_unit' must be provided together".to_string(),
                ));
            }
        };

        Ok(Self {
            factor,
            offset,
            unit_conversion,
        })
    }
}

/// Meters per one linear unit for the supported presets.
fn meters_per_unit(unit: &str) -> Result<f64, ToolError> {
    let u = unit.trim().to_ascii_lowercase();
    let m = match u.as_str() {
        "meters" | "meter" | "metre" | "m" => 1.0,
        "feet" | "foot" | "ft" => 0.3048, // international foot
        "us_feet" | "us_foot" | "us_survey_foot" | "us_survey_feet" => 1200.0 / 3937.0,
        "centimeters" | "centimeter" | "cm" => 0.01,
        "millimeters" | "millimeter" | "mm" => 0.001,
        "kilometers" | "kilometer" | "km" => 1000.0,
        "miles" | "mile" | "mi" => 1609.344,
        "yards" | "yard" | "yd" => 0.9144,
        "inches" | "inch" | "in" => 0.0254,
        other => {
            return Err(ToolError::Validation(format!(
                "unknown unit '{other}' (expected one of: meters, feet, us_feet, centimeters, millimeters, kilometers, miles, yards, inches)"
            )));
        }
    };
    Ok(m)
}

/// Parses an optional numeric parameter accepting a JSON number OR a numeric
/// string (host UIs post scalars as strings). Absent / null / empty -> None.
fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

/// Applies `f` to every coordinate in `geom`, mutating it in place.
fn for_each_coord_mut(geom: &mut Geometry, f: &mut impl FnMut(&mut Coord)) {
    match geom {
        Geometry::Point(c) => f(c),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => {
            for c in cs.iter_mut() {
                f(c);
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            for c in exterior.0.iter_mut() {
                f(c);
            }
            for r in interiors.iter_mut() {
                for c in r.0.iter_mut() {
                    f(c);
                }
            }
        }
        Geometry::MultiLineString(ls) => {
            for c in ls.iter_mut().flatten() {
                f(c);
            }
        }
        Geometry::MultiPolygon(ps) => {
            for (e, hs) in ps.iter_mut() {
                for c in e.0.iter_mut() {
                    f(c);
                }
                for r in hs.iter_mut() {
                    for c in r.0.iter_mut() {
                        f(c);
                    }
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs.iter_mut() {
                for_each_coord_mut(g, f);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Layer, Ring};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = Adjust3dZTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn z_layer() -> String {
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(
                Some(Geometry::point_z(0.0, 0.0, 100.0)),
                &[("name", "a".into())],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::point_z(1.0, 1.0, 200.0)),
                &[("name", "b".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn first_z(geom: &Geometry) -> Option<f64> {
        geom.all_coords().first().and_then(|c| c.z)
    }

    #[test]
    fn feet_to_meters_scales_z() {
        let input = z_layer();
        let (out, layer) = run_tool(json!({
            "input": input,
            "from_unit": "feet",
            "to_unit": "meters",
        }));
        assert_eq!(out.outputs["features_with_z"], json!(2));
        assert_eq!(out.outputs["vertices_adjusted"], json!(2));
        // 100 ft -> 30.48 m, 200 ft -> 60.96 m
        let z0 = first_z(layer.features[0].geometry.as_ref().unwrap()).unwrap();
        assert!((z0 - 30.48).abs() < 1e-9, "z0 = {z0}");
        let z1 = first_z(layer.features[1].geometry.as_ref().unwrap()).unwrap();
        assert!((z1 - 60.96).abs() < 1e-9, "z1 = {z1}");
        assert!((out.outputs["z_max"].as_f64().unwrap() - 60.96).abs() < 1e-9);
    }

    #[test]
    fn factor_and_offset_apply() {
        let input = z_layer();
        let (_out, layer) = run_tool(json!({
            "input": input,
            "factor": 2.0,
            "offset": 5.0,
        }));
        // 100 * 2 + 5 = 205
        let z0 = first_z(layer.features[0].geometry.as_ref().unwrap()).unwrap();
        assert!((z0 - 205.0).abs() < 1e-9, "z0 = {z0}");
    }

    #[test]
    fn string_params_accepted() {
        let input = z_layer();
        let (_out, layer) = run_tool(json!({
            "input": input,
            "factor": "0.5",
            "offset": "10",
        }));
        // 100 * 0.5 + 10 = 60
        let z0 = first_z(layer.features[0].geometry.as_ref().unwrap()).unwrap();
        assert!((z0 - 60.0).abs() < 1e-9, "z0 = {z0}");
    }

    #[test]
    fn polygon_z_all_vertices_adjusted() {
        let mut layer = Layer::new("poly");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        let ext = Ring::new(vec![
            Coord::xyz(0.0, 0.0, 10.0),
            Coord::xyz(4.0, 0.0, 10.0),
            Coord::xyz(4.0, 4.0, 10.0),
            Coord::xyz(0.0, 4.0, 10.0),
        ]);
        layer
            .add_feature(
                Some(Geometry::Polygon {
                    exterior: ext,
                    interiors: vec![],
                }),
                &[("name", "sq".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "offset": 90.0 }));
        assert_eq!(out.outputs["vertices_adjusted"], json!(4));
        for c in layer.features[0].geometry.as_ref().unwrap().all_coords() {
            assert!((c.z.unwrap() - 100.0).abs() < 1e-9);
        }
    }

    #[test]
    fn features_without_z_pass_through() {
        let mut layer = Layer::new("flat");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(Some(Geometry::point(3.0, 4.0)), &[("name", "flat".into())])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "factor": 3.0, "offset": 1.0 }));
        assert_eq!(out.outputs["features_with_z"], json!(0));
        assert_eq!(out.outputs["vertices_adjusted"], json!(0));
        // 2D geometry unchanged, no z introduced.
        let c = layer.features[0].geometry.as_ref().unwrap().all_coords()[0].clone();
        assert_eq!((c.x, c.y), (3.0, 4.0));
        assert!(c.z.is_none());
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = z_layer();
        // factor not a number
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input.clone(),
            "factor": "abc",
        }))
        .unwrap();
        assert!(Adjust3dZTool.run(&args, &ctx()).is_err());

        // only one of from_unit/to_unit
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input.clone(),
            "from_unit": "feet",
        }))
        .unwrap();
        assert!(Adjust3dZTool.run(&args, &ctx()).is_err());

        // unknown unit
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input,
            "from_unit": "furlongs",
            "to_unit": "meters",
        }))
        .unwrap();
        assert!(Adjust3dZTool.run(&args, &ctx()).is_err());
    }
}
