//! GeoLibre tool: volumetric union of overlapping 3D solids.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Union 3D* (3D Analyst).
//!
//! The 2D overlay suite is complete and battle-tested — `union`, `intersect`,
//! `erase`, `identity`, `symmetrical_difference`, `clip`, all backed by `geo`
//! `BooleanOps`. There is no 3D counterpart anywhere in either registry.
//!
//! That gap has a concrete cost: the moment solids overlap, **volumes stop being
//! additive**. Summing `minimum_bounding_volume` outputs, or `buffer_3d`
//! capsules around a pipe network, or per-source plume envelopes, double-counts
//! every overlap. `polygon_volume` and `surface_volume` both measure against a
//! *reference plane* and cannot help with solid-solid overlap. So there was no
//! correct way to answer "what total volume do these solids cover".
//!
//! ## Scope, deliberately
//!
//! This ships the **volume answer**, not an exact mesh boolean. Per-solid
//! volumes are exact (signed-tetrahedron summation). The *union* volume is
//! computed by voxel occupancy over a shared grid: approximate, but bounded,
//! trivially parallel, and free of the floating-point robustness minefield that
//! an exact mesh boolean's arrangement/coplanar-face handling represents. Most
//! callers want the number rather than the merged mesh, and the accuracy is
//! reported (`resolution`, plus a convergence-friendly knob) instead of being
//! implied.
//!
//! Exact merged geometry is intentionally left for a follow-up; the parameter
//! surface here does not promise it.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::inside_3d::{collect_triangles, Solid, Tri};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Computes overlap-corrected combined volumes for groups of 3D solids.
pub struct Union3dTool;

impl Tool for Union3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "union_3d",
            display_name: "Union 3D",
            summary: "Computes the combined volume of overlapping closed 3D solids without double-counting their intersections, plus per-pair overlap volumes (ArcGIS Union 3D). Neither registry has any 3D overlay, so summing buffer_3d capsules or minimum_bounding_volume hulls today double-counts every overlap; polygon_volume and surface_volume only measure against a reference plane. Per-solid volumes are exact (signed-tetrahedron summation); the union is estimated by voxel occupancy at a reported resolution.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Closed 3D solid features (triangle-mesh MultiPolygons with Z, as buffer_3d and minimum_bounding_volume emit).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output table path of per-group volumes. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_field",
                    description: "Optional field; solids are unioned within each group value rather than across the whole layer.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_overlap_table",
                    description: "Optional path for the per-pair overlap volumes. If omitted, stored in memory (still returned).",
                    required: false,
                },
                ToolParamSpec {
                    name: "resolution",
                    description: "Voxel samples per axis used to estimate the union volume, 4 to 512 (default 96). Higher is more accurate and slower.",
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
        parse_resolution(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = required_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let overlap_out = parse_optional_str(args, "output_overlap_table")?;
        let group_field = parse_optional_str(args, "group_field")?;
        let resolution = parse_resolution(args)?;

        let layer = load_input_layer(input)?;
        let group_idx = match group_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("group_field '{f}' not found on the input layer"))
            })?),
            None => None,
        };

        // Load solids, grouped.
        let mut groups: BTreeMap<String, Vec<Solid>> = BTreeMap::new();
        let mut open_meshes = 0_u64;
        for (fid, feature) in layer.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let tris = collect_triangles(geom);
            if tris.is_empty() {
                continue;
            }
            let solid = Solid::new(fid, tris);
            if !solid.closed {
                // An open mesh does not bound a volume; its "volume" would be
                // meaningless, so report and skip rather than emit a number.
                open_meshes += 1;
                continue;
            }
            let key = match group_idx {
                Some(i) => field_to_string(&feature.attributes[i]),
                None => "all".to_string(),
            };
            groups.entry(key).or_default().push(solid);
        }
        if groups.is_empty() {
            return Err(ToolError::Execution(format!(
                "input holds no closed triangle-mesh solids ({open_meshes} open mesh(es) skipped)"
            )));
        }
        ctx.progress
            .info(&format!("{} group(s) of solids", groups.len()));

        let mut out = Layer::new("union_3d");
        out.add_field(FieldDef::new("group", FieldType::Text));
        out.add_field(FieldDef::new("solid_count", FieldType::Integer));
        out.add_field(FieldDef::new("sum_volume", FieldType::Float));
        out.add_field(FieldDef::new("union_volume", FieldType::Float));
        out.add_field(FieldDef::new("overlap_volume", FieldType::Float));

        let mut overlap = Layer::new("union_3d_overlaps");
        overlap.add_field(FieldDef::new("group", FieldType::Text));
        overlap.add_field(FieldDef::new("fid_a", FieldType::Integer));
        overlap.add_field(FieldDef::new("fid_b", FieldType::Integer));
        overlap.add_field(FieldDef::new("overlap_volume", FieldType::Float));

        let mut total_union = 0.0_f64;
        let mut total_sum = 0.0_f64;
        let mut overlap_pairs = 0_u64;
        let mut total_sampled = 0_u64;

        for (gi, (key, solids)) in groups.iter().enumerate() {
            // Exact per-solid volumes.
            let vols: Vec<f64> = solids.iter().map(|s| mesh_volume(&s.tris)).collect();
            let sum_volume: f64 = vols.iter().sum();

            // Split into bbox-connected components: an isolated solid keeps its
            // exact volume, and only genuinely interacting clusters pay for
            // sampling. Without this a single touching pair would downgrade an
            // entire group of otherwise disjoint solids to an estimate, and
            // over a needlessly large grid.
            let mut union_volume = 0.0_f64;
            let mut sampled_components = 0_u64;
            for comp in bbox_components(solids) {
                if comp.len() == 1 {
                    union_volume += vols[comp[0]];
                } else {
                    let members: Vec<&Solid> = comp.iter().map(|i| &solids[*i]).collect();
                    union_volume += occupancy_volume(&members, resolution);
                    sampled_components += 1;
                }
            }

            total_sum += sum_volume;
            total_union += union_volume;
            total_sampled += sampled_components;

            out.add_feature(
                None,
                &[
                    ("group", FieldValue::Text(key.clone())),
                    ("solid_count", FieldValue::Integer(solids.len() as i64)),
                    ("sum_volume", FieldValue::Float(sum_volume)),
                    ("union_volume", FieldValue::Float(union_volume)),
                    (
                        "overlap_volume",
                        // What naive summation would have double-counted.
                        FieldValue::Float((sum_volume - union_volume).max(0.0)),
                    ),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing group row: {e}")))?;

            // Pairwise overlaps, for the pairs whose boxes actually meet.
            for a in 0..solids.len() {
                for b in (a + 1)..solids.len() {
                    if !bbox_overlap(&solids[a], &solids[b]) {
                        continue;
                    }
                    let v = pair_intersection_volume(&solids[a], &solids[b], resolution);
                    if v <= 0.0 {
                        continue;
                    }
                    overlap_pairs += 1;
                    overlap
                        .add_feature(
                            None,
                            &[
                                ("group", FieldValue::Text(key.clone())),
                                ("fid_a", FieldValue::Integer(solids[a].fid as i64)),
                                ("fid_b", FieldValue::Integer(solids[b].fid as i64)),
                                ("overlap_volume", FieldValue::Float(v)),
                            ],
                        )
                        .map_err(|e| {
                            ToolError::Execution(format!("failed writing overlap row: {e}"))
                        })?;
                }
            }

            ctx.progress
                .progress((gi as f64 + 1.0) / groups.len() as f64);
        }

        let group_count = groups.len();
        let out_path = write_or_store_layer(out, output)?;
        let overlap_path = write_or_store_layer(overlap, overlap_out)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_overlap_table".to_string(), json!(overlap_path));
        outputs.insert("group_count".to_string(), json!(group_count));
        outputs.insert("sum_volume".to_string(), json!(total_sum));
        outputs.insert("union_volume".to_string(), json!(total_union));
        outputs.insert(
            "overlap_volume".to_string(),
            json!((total_sum - total_union).max(0.0)),
        );
        outputs.insert("open_mesh_count".to_string(), json!(open_meshes));
        outputs.insert("overlap_pair_count".to_string(), json!(overlap_pairs));
        outputs.insert("resolution".to_string(), json!(resolution));
        outputs.insert("sampled_component_count".to_string(), json!(total_sampled));
        Ok(ToolRunResult { outputs })
    }
}

/// Exact volume of a closed triangle mesh by signed-tetrahedron summation.
fn mesh_volume(tris: &[Tri]) -> f64 {
    let mut v = 0.0;
    for t in tris {
        let (a, b, c) = (t[0], t[1], t[2]);
        v += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    (v / 6.0).abs()
}

fn bbox_overlap(a: &Solid, b: &Solid) -> bool {
    (0..3).all(|k| a.min[k] <= b.max[k] && a.max[k] >= b.min[k])
}

/// Groups solid indices into components whose bounding boxes transitively
/// overlap. Singletons are exactly measurable; only multi-member components
/// need sampling.
fn bbox_components(solids: &[Solid]) -> Vec<Vec<usize>> {
    let n = solids.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if bbox_overlap(&solids[i], &solids[j]) {
                let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                if a != b {
                    parent[a] = b;
                }
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    groups.into_values().collect()
}

/// Union volume by voxel occupancy: a cell counts once if it is inside **any**
/// solid, which is exactly what stops overlaps being double-counted.
fn occupancy_volume(solids: &[&Solid], resolution: usize) -> f64 {
    let Some((min, max)) = union_bbox(solids) else {
        return 0.0;
    };
    let (nx, ny, nz, cell_vol, step) = grid_for(min, max, resolution);
    if cell_vol <= 0.0 {
        return 0.0;
    }
    let mut occupied = 0_u64;
    for k in 0..nz {
        let z = min[2] + (k as f64 + 0.5) * step[2];
        for j in 0..ny {
            let y = min[1] + (j as f64 + 0.5) * step[1];
            for i in 0..nx {
                let x = min[0] + (i as f64 + 0.5) * step[0];
                if solids.iter().any(|s| s.contains(x, y, z)) {
                    occupied += 1;
                }
            }
        }
    }
    occupied as f64 * cell_vol
}

/// Intersection volume of one pair, sampled over their shared bounding box.
fn pair_intersection_volume(a: &Solid, b: &Solid, resolution: usize) -> f64 {
    let mut min = [0.0_f64; 3];
    let mut max = [0.0_f64; 3];
    for k in 0..3 {
        min[k] = a.min[k].max(b.min[k]);
        max[k] = a.max[k].min(b.max[k]);
        if max[k] <= min[k] {
            return 0.0;
        }
    }
    let (nx, ny, nz, cell_vol, step) = grid_for(min, max, resolution);
    if cell_vol <= 0.0 {
        return 0.0;
    }
    let mut both = 0_u64;
    for k in 0..nz {
        let z = min[2] + (k as f64 + 0.5) * step[2];
        for j in 0..ny {
            let y = min[1] + (j as f64 + 0.5) * step[1];
            for i in 0..nx {
                let x = min[0] + (i as f64 + 0.5) * step[0];
                if a.contains(x, y, z) && b.contains(x, y, z) {
                    both += 1;
                }
            }
        }
    }
    both as f64 * cell_vol
}

fn union_bbox(solids: &[&Solid]) -> Option<([f64; 3], [f64; 3])> {
    if solids.is_empty() {
        return None;
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for s in solids {
        for k in 0..3 {
            min[k] = min[k].min(s.min[k]);
            max[k] = max[k].max(s.max[k]);
        }
    }
    if (0..3).any(|k| max[k] <= min[k]) {
        return None;
    }
    Some((min, max))
}

/// Builds a sampling grid over a box: the longest axis gets `resolution` cells
/// and the others are scaled to keep voxels cubic, so accuracy does not depend
/// on the box's aspect ratio.
fn grid_for(
    min: [f64; 3],
    max: [f64; 3],
    resolution: usize,
) -> (usize, usize, usize, f64, [f64; 3]) {
    let span = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
    let longest = span[0].max(span[1]).max(span[2]);
    if longest <= 0.0 {
        return (0, 0, 0, 0.0, [0.0; 3]);
    }
    let h = longest / resolution as f64;
    let n: Vec<usize> = span
        .iter()
        .map(|s| ((s / h).ceil() as usize).max(1))
        .collect();
    let step = [
        span[0] / n[0] as f64,
        span[1] / n[1] as f64,
        span[2] / n[2] as f64,
    ];
    (n[0], n[1], n[2], step[0] * step[1] * step[2], step)
}

fn field_to_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null => String::new(),
        other => format!("{other:?}"),
    }
}

fn parse_resolution(args: &ToolArgs) -> Result<usize, ToolError> {
    let v = match args.get("resolution") {
        None | Some(Value::Null) => return Ok(96),
        Some(Value::Number(n)) => n.as_f64().unwrap_or(f64::NAN),
        Some(Value::String(s)) if s.trim().is_empty() => return Ok(96),
        Some(Value::String(s)) => s.trim().parse::<f64>().map_err(|_| {
            ToolError::Validation("parameter 'resolution' must be a number".to_string())
        })?,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'resolution' must be a number".to_string(),
            ))
        }
    };
    if !v.is_finite() || !(4.0..=512.0).contains(&v) {
        return Err(ToolError::Validation(
            "'resolution' must be between 4 and 512 samples per axis".to_string(),
        ));
    }
    Ok(v as usize)
}

fn required_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::inside_3d::box_mesh;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn solids(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("solids");
        l.geom_type = Some(GeometryType::MultiPolygon);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(extra: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(extra).unwrap();
        let res = Union3dTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (l, res)
    }

    /// A single cube's volume is exact — signed-tetrahedron summation, not a
    /// sampled estimate.
    #[test]
    fn single_solid_volume_is_exact() {
        let path = solids(vec![box_mesh([0.0, 0.0, 0.0], [2.0, 3.0, 4.0])]);
        let (_, res) = run(json!({ "input": path }));
        let v = res.outputs["union_volume"].as_f64().unwrap();
        assert!((v - 24.0).abs() < 1e-9, "expected exactly 24, got {v}");
        assert_eq!(res.outputs["overlap_volume"].as_f64().unwrap(), 0.0);
    }

    /// Disjoint solids are additive, and take the exact path rather than the
    /// sampled one.
    #[test]
    fn disjoint_solids_are_additive_and_exact() {
        let path = solids(vec![
            box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_mesh([10.0, 0.0, 0.0], [12.0, 1.0, 1.0]),
        ]);
        let (_, res) = run(json!({ "input": path }));
        let v = res.outputs["union_volume"].as_f64().unwrap();
        assert!((v - 3.0).abs() < 1e-9, "expected exactly 3, got {v}");
        assert_eq!(res.outputs["overlap_volume"].as_f64().unwrap(), 0.0);
    }

    /// The whole reason the tool exists: overlapping solids must NOT be summed.
    /// Two unit cubes overlapping in half their volume cover 1.5, not 2.
    #[test]
    fn overlapping_solids_are_not_double_counted() {
        let path = solids(vec![
            box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_mesh([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ]);
        let (_, res) = run(json!({ "input": path, "resolution": 120 }));
        let sum = res.outputs["sum_volume"].as_f64().unwrap();
        let union = res.outputs["union_volume"].as_f64().unwrap();
        assert!((sum - 2.0).abs() < 1e-9, "naive sum is 2, got {sum}");
        assert!(
            (union - 1.5).abs() < 0.02,
            "union of the two half-overlapping cubes is 1.5, got {union}"
        );
        let ov = res.outputs["overlap_volume"].as_f64().unwrap();
        assert!((ov - 0.5).abs() < 0.02, "overlap is 0.5, got {ov}");
    }

    /// Fully nested solids: the union is just the outer volume.
    #[test]
    fn nested_solid_contributes_nothing() {
        let path = solids(vec![
            box_mesh([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]),
            box_mesh([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
        ]);
        let (_, res) = run(json!({ "input": path, "resolution": 100 }));
        let union = res.outputs["union_volume"].as_f64().unwrap();
        assert!(
            (union - 64.0).abs() < 1.5,
            "a fully nested cube adds nothing; expected ~64, got {union}"
        );
    }

    /// The pairwise overlap table records the shared volume.
    #[test]
    fn pairwise_overlap_is_reported() {
        let path = solids(vec![
            box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_mesh([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "resolution": 120 })).unwrap();
        let res = Union3dTool.run(&args, &ctx()).unwrap();
        assert_eq!(res.outputs["overlap_pair_count"], json!(1));

        let table =
            load_input_layer(res.outputs["output_overlap_table"].as_str().unwrap()).unwrap();
        assert_eq!(table.len(), 1);
        let i = table.schema.field_index("overlap_volume").unwrap();
        let FieldValue::Float(v) = table.iter().next().unwrap().attributes[i] else {
            panic!("expected a float overlap volume");
        };
        assert!((v - 0.5).abs() < 0.02, "pair overlap is 0.5, got {v}");
    }

    /// A single overlapping pair must not downgrade unrelated disjoint solids in
    /// the same group to a sampled estimate.
    #[test]
    fn disjoint_solids_stay_exact_alongside_an_overlapping_pair() {
        let path = solids(vec![
            // One overlapping pair...
            box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_mesh([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]),
            // ...plus two solids far away from everything.
            box_mesh([100.0, 0.0, 0.0], [102.0, 1.0, 1.0]),
            box_mesh([200.0, 0.0, 0.0], [203.0, 1.0, 1.0]),
        ]);
        let (_, res) = run(json!({ "input": path, "resolution": 64 }));
        // Only the touching pair needs sampling.
        assert_eq!(res.outputs["sampled_component_count"], json!(1));
        // 1.75 (pair) + 2 + 3, with the two isolated boxes contributing exactly.
        let v = res.outputs["union_volume"].as_f64().unwrap();
        assert!((v - 6.75).abs() < 0.05, "expected ~6.75, got {v}");
    }

    /// Grouping unions within each group value independently.
    #[test]
    fn group_field_partitions_the_union() {
        let mut l = Layer::new("solids");
        l.geom_type = Some(GeometryType::MultiPolygon);
        l.add_field(FieldDef::new("site", FieldType::Text));
        // Two overlapping cubes on site A, one separate cube on site B.
        l.add_feature(
            Some(box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])),
            &[("site", FieldValue::Text("A".into()))],
        )
        .unwrap();
        l.add_feature(
            Some(box_mesh([0.5, 0.0, 0.0], [1.5, 1.0, 1.0])),
            &[("site", FieldValue::Text("A".into()))],
        )
        .unwrap();
        l.add_feature(
            Some(box_mesh([50.0, 0.0, 0.0], [52.0, 1.0, 1.0])),
            &[("site", FieldValue::Text("B".into()))],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);

        let (table, res) = run(json!({
            "input": path, "group_field": "site", "resolution": 120
        }));
        assert_eq!(res.outputs["group_count"], json!(2));
        assert_eq!(table.len(), 2);

        let gi = table.schema.field_index("group").unwrap();
        let vi = table.schema.field_index("union_volume").unwrap();
        for f in table.iter() {
            let FieldValue::Text(g) = &f.attributes[gi] else {
                panic!("expected text group");
            };
            let FieldValue::Float(v) = f.attributes[vi] else {
                panic!("expected float volume");
            };
            match g.as_str() {
                "A" => assert!((v - 1.5).abs() < 0.03, "site A union is 1.5, got {v}"),
                "B" => assert!(
                    (v - 2.0).abs() < 1e-9,
                    "site B is a lone exact cube, got {v}"
                ),
                other => panic!("unexpected group {other}"),
            }
        }
    }

    /// Raising the resolution moves the estimate toward the true value, which is
    /// what makes the approximation honest rather than arbitrary.
    ///
    /// The two boxes are offset on **two** axes so their union is an L-shaped
    /// prism. Offsetting on a single axis would make the union itself a box,
    /// which voxel occupancy reproduces exactly at any resolution — the
    /// comparison would then say nothing.
    #[test]
    fn higher_resolution_converges() {
        let path = solids(vec![
            box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_mesh([0.5, 0.5, 0.0], [1.5, 1.5, 1.0]),
        ]);
        // 1 + 1 - (0.5 * 0.5 * 1) = 1.75
        let truth = 1.75;
        let (_, coarse) = run(json!({ "input": path.clone(), "resolution": 5 }));
        let (_, fine) = run(json!({ "input": path, "resolution": 200 }));
        let e_coarse = (coarse.outputs["union_volume"].as_f64().unwrap() - truth).abs();
        let e_fine = (fine.outputs["union_volume"].as_f64().unwrap() - truth).abs();
        assert!(
            e_fine < e_coarse,
            "finer sampling must be closer: coarse err {e_coarse}, fine err {e_fine}"
        );
        assert!(
            e_fine < 0.01,
            "at resolution 200 the estimate should be within 1%, err {e_fine}"
        );
    }

    /// Open meshes bound no volume and are skipped, not silently measured.
    #[test]
    fn open_meshes_are_skipped_and_counted() {
        let open = Geometry::MultiPolygon(vec![(
            wbvector::Ring::new(vec![
                wbvector::Coord::xyz(0.0, 0.0, 0.0),
                wbvector::Coord::xyz(1.0, 0.0, 0.0),
                wbvector::Coord::xyz(0.0, 1.0, 0.0),
            ]),
            Vec::new(),
        )]);
        let path = solids(vec![box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), open]);
        let (_, res) = run(json!({ "input": path }));
        assert_eq!(res.outputs["open_mesh_count"], json!(1));
        let v = res.outputs["union_volume"].as_f64().unwrap();
        assert!(
            (v - 1.0).abs() < 1e-9,
            "only the closed cube counts, got {v}"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(Union3dTool.validate(&args).is_err());

        let path = solids(vec![box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])]);
        // A blank input must fail validation, not only at run time.
        let blank: ToolArgs = serde_json::from_value(json!({ "input": "   " })).unwrap();
        assert!(Union3dTool.validate(&blank).is_err());

        for bad in [
            json!({ "input": path.clone(), "resolution": 2 }),
            json!({ "input": path.clone(), "resolution": 5000 }),
            json!({ "input": path.clone(), "resolution": "lots" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(Union3dTool.validate(&args).is_err());
        }

        // A missing group_field is a run-time error.
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "group_field": "nope" })).unwrap();
        assert!(Union3dTool.run(&args, &ctx()).is_err());
    }
}
