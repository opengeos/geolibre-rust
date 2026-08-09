//! GeoLibre tool: build a spectral library from an image plus training regions.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Extract Spectra From Image*
//! (Image Analyst).
//!
//! ## The producer the spectral suite never had
//!
//! Three shipped tools **consume** a spectral library and nothing created one:
//! `spectral_library_matching`, `spectral_angle_mapper`, and
//! `matched_filter_target_detection` (`linear_spectral_unmixing` likewise needs
//! endmember spectra). A user therefore had to hand-author the library file,
//! which in practice meant the whole spectral suite only worked against a
//! pre-existing published library, never against their own scene.
//!
//! ## Output format is the requirement, not a detail
//!
//! The library is written in exactly the layout the bundled
//! `spectral_library_matching` parses for its `library_csv` parameter:
//! one headerless line per class, `name,b1,b2,...`, with **exactly one value per
//! band** and `#` comments permitted. A library this tool emits that its
//! consumer cannot read would make the whole thing pointless, so the per-class
//! standard deviations and cell counts go to a **separate** `output_stats` CSV
//! rather than widening the library rows — extra columns there are a hard parse
//! error downstream.
//!
//! ## Why trimming matters
//!
//! Training polygons drawn by hand always clip a few mixed edge pixels, and a
//! single bright roof inside a "vegetation" polygon shifts the mean endmember
//! enough to mis-classify a scene. `trim_percent` discards that fraction from
//! each tail per band before reducing, which is the cheap standard defence.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{FieldValue, Geometry};

use crate::args_common::{choice_or, opt_f64, req_str, usize_or};
use crate::common::{load_input_raster, write_text_output};
use crate::sar_common::envelope;
use crate::vector_common::{geometry_contains_point, load_input_layer, parse_optional_str};

const STATISTICS: [&str; 4] = ["mean", "median", "max", "min"];

pub struct ExtractSpectraFromImageTool;

impl Tool for ExtractSpectraFromImageTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "extract_spectra_from_image",
            display_name: "Extract Spectra From Image",
            summary: "Derives a spectral library from a multiband image and training regions, one spectrum per class, in the CSV layout spectral_library_matching reads (ArcGIS Extract Spectra From Image). spectral_library_matching, spectral_angle_mapper and matched_filter_target_detection all consume a library and nothing created one, so the spectral suite could only be used against pre-existing published libraries.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Multiband image.",
                    required: true,
                },
                ToolParamSpec {
                    name: "training_features",
                    description: "Polygon or point layer delineating each class or target.",
                    required: true,
                },
                ToolParamSpec {
                    name: "class_field",
                    description: "Attribute holding the class label. Features sharing a label are pooled into one spectrum.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output spectral library CSV path ('name,b1,b2,...' per class), ready for spectral_library_matching's 'library_csv'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "Reduction across the covered cells: 'mean' (default), 'median', 'max', 'min'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "trim_percent",
                    description: "Percentage discarded from each tail per band before reducing, to suppress mixed edge pixels (default 0, max 49).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_cells",
                    description: "Classes covering fewer cells than this are reported and skipped rather than emitting a one-pixel 'spectrum' (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "wavelengths",
                    description: "Optional comma-separated per-band wavelengths, recorded as a comment header in the library and in output_stats.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_stats",
                    description: "Optional CSV receiving per-class, per-band standard deviation and cell count. Kept separate so the library stays parseable.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "training_features")?;
        req_str(args, "class_field")?;
        req_str(args, "output")?;
        choice_or(args, "statistic", &STATISTICS, "mean")?;
        if let Some(t) = opt_f64(args, "trim_percent")? {
            // 50% from each tail would discard every value.
            if !(0.0..=49.0).contains(&t) {
                return Err(ToolError::Validation(
                    "'trim_percent' must be between 0 and 49 (it is applied to EACH tail)"
                        .to_string(),
                ));
            }
        }
        if usize_or(args, "min_cells", 1)? == 0 {
            return Err(ToolError::Validation(
                "'min_cells' must be at least 1".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let training = req_str(args, "training_features")?.to_string();
        let class_field = req_str(args, "class_field")?.to_string();
        let output = req_str(args, "output")?.to_string();
        let statistic = choice_or(args, "statistic", &STATISTICS, "mean")?;
        let trim = opt_f64(args, "trim_percent")?.unwrap_or(0.0);
        let min_cells = usize_or(args, "min_cells", 1)?.max(1);
        let output_stats = parse_optional_str(args, "output_stats")?.map(str::to_string);

        let raster = load_input_raster(&input)?;
        let bands = raster.bands;
        if bands == 0 {
            return Err(ToolError::Execution(
                "input raster has no bands".to_string(),
            ));
        }
        let wavelengths = parse_wavelengths(args, "wavelengths", bands)?;

        let layer = load_input_layer(&training)?;
        let fidx = layer.schema.field_index(&class_field).ok_or_else(|| {
            ToolError::Validation(format!(
                "class_field '{class_field}' not found in the training layer"
            ))
        })?;

        // Per class, per band, the raw covered-cell values. Held rather than
        // accumulated because median and trimming both need the full vector.
        let mut pools: BTreeMap<String, Vec<Vec<f64>>> = BTreeMap::new();
        let mut cell_counts: BTreeMap<String, usize> = BTreeMap::new();

        let rows = raster.rows;
        let cols = raster.cols;
        let y_max = raster.y_min + rows as f64 * raster.cell_size_y;

        for feature in layer.iter() {
            let Some(label) = feature.attributes.get(fidx).map(label_string) else {
                continue;
            };
            if label.is_empty() {
                continue;
            }
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };

            let entry = pools
                .entry(label.clone())
                .or_insert_with(|| vec![Vec::new(); bands]);
            let count = cell_counts.entry(label).or_insert(0);

            match geom {
                // Points sample the single cell they fall in — nearest cell, no
                // interpolation, because a spectrum is a per-pixel measurement.
                Geometry::Point(p) => {
                    if let Some((r, c)) = cell_of(&raster, y_max, p.x, p.y) {
                        if push_cell(&raster, bands, r, c, entry) {
                            *count += 1;
                        }
                    }
                }
                Geometry::MultiPoint(ps) => {
                    for p in ps {
                        if let Some((r, c)) = cell_of(&raster, y_max, p.x, p.y) {
                            if push_cell(&raster, bands, r, c, entry) {
                                *count += 1;
                            }
                        }
                    }
                }
                _ => {
                    // Bounding-box prefilter, then the exact predicate. The
                    // prefilter must cover exactly the variants the predicate
                    // covers (round-18 lesson: sar_common::envelope handles
                    // GeometryCollection, and geometry_contains_point recurses
                    // into it, so the two agree).
                    // An empty geometry yields an inverted box, which makes
                    // the loop bounds below empty — the same answer the
                    // containment test would give, only far more cheaply.
                    let (min_x, min_y, max_x, max_y) = envelope(geom);
                    let c0 = (((min_x - raster.x_min) / raster.cell_size_x).floor() as isize)
                        .clamp(0, cols as isize - 1) as usize;
                    let c1 = (((max_x - raster.x_min) / raster.cell_size_x).ceil() as isize)
                        .clamp(0, cols as isize - 1) as usize;
                    let r0 = (((y_max - max_y) / raster.cell_size_y).floor() as isize)
                        .clamp(0, rows as isize - 1) as usize;
                    let r1 = (((y_max - min_y) / raster.cell_size_y).ceil() as isize)
                        .clamp(0, rows as isize - 1) as usize;
                    for r in r0..=r1 {
                        for c in c0..=c1 {
                            let x = raster.x_min + (c as f64 + 0.5) * raster.cell_size_x;
                            let y = y_max - (r as f64 + 0.5) * raster.cell_size_y;
                            if geometry_contains_point(geom, x, y)
                                && push_cell(&raster, bands, r, c, entry)
                            {
                                *count += 1;
                            }
                        }
                    }
                }
            }
        }

        if pools.is_empty() {
            return Err(ToolError::Execution(
                "no training feature covered a valid image cell".to_string(),
            ));
        }

        let mut lines: Vec<String> = Vec::new();
        if let Some(w) = &wavelengths {
            // A '#' comment: the parser skips these, so it is safe metadata.
            lines.push(format!(
                "# wavelengths: {}",
                w.iter().map(fmt).collect::<Vec<_>>().join(",")
            ));
        }
        let mut stats_lines: Vec<String> = vec!["class,band,wavelength,value,stddev,count".into()];
        let mut skipped: Vec<String> = Vec::new();
        let mut emitted: Vec<String> = Vec::new();

        for (label, per_band) in &pools {
            let n = cell_counts.get(label).copied().unwrap_or(0);
            if n < min_cells {
                skipped.push(label.clone());
                continue;
            }
            let mut values = Vec::with_capacity(bands);
            for (b, raw) in per_band.iter().enumerate() {
                let mut v = raw.clone();
                v.sort_by(f64::total_cmp);
                let v = trim_tails(&v, trim);
                let value = reduce(&v, statistic);
                let sd = stddev(&v);
                values.push(value);
                stats_lines.push(format!(
                    "{},{},{},{},{},{}",
                    csv_escape(label),
                    b + 1,
                    wavelengths.as_ref().map(|w| fmt(&w[b])).unwrap_or_default(),
                    fmt(&value),
                    fmt(&sd),
                    v.len()
                ));
            }
            lines.push(format!(
                "{},{}",
                csv_escape(label),
                values.iter().map(fmt).collect::<Vec<_>>().join(",")
            ));
            emitted.push(label.clone());
        }

        if emitted.is_empty() {
            return Err(ToolError::Execution(format!(
                "every class covered fewer than min_cells ({min_cells}) valid cells; lower \
                 'min_cells' or enlarge the training regions"
            )));
        }
        ctx.progress.info(&format!(
            "{} class(es) written, {} skipped, {bands} band(s)",
            emitted.len(),
            skipped.len()
        ));

        write_text_output(&format!("{}\n", lines.join("\n")), &output)?;
        let stats_path = match &output_stats {
            Some(p) => {
                write_text_output(&format!("{}\n", stats_lines.join("\n")), p)?;
                Some(p.clone())
            }
            None => None,
        };

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        if let Some(p) = stats_path {
            outputs.insert("output_stats".to_string(), json!(p));
        }
        outputs.insert("class_count".to_string(), json!(emitted.len()));
        outputs.insert("band_count".to_string(), json!(bands));
        outputs.insert("classes".to_string(), json!(emitted));
        outputs.insert("skipped_classes".to_string(), json!(skipped));
        outputs.insert("statistic".to_string(), json!(statistic));
        Ok(ToolRunResult { outputs })
    }
}

/// Appends a cell's per-band values to the pool.
///
/// Returns false (contributing nothing) when **any** band is no-data at that
/// cell: a pixel missing one band is not a valid spectrum, and letting it
/// through would make the bands of one "spectrum" come from different cells.
fn push_cell(raster: &Raster, bands: usize, r: usize, c: usize, entry: &mut [Vec<f64>]) -> bool {
    let mut vals = Vec::with_capacity(bands);
    for b in 0..bands {
        let v = raster.get(b as isize, r as isize, c as isize);
        if v == raster.nodata || !v.is_finite() {
            return false;
        }
        vals.push(v);
    }
    for (b, v) in vals.into_iter().enumerate() {
        entry[b].push(v);
    }
    true
}

fn cell_of(raster: &Raster, y_max: f64, x: f64, y: f64) -> Option<(usize, usize)> {
    let c = ((x - raster.x_min) / raster.cell_size_x).floor();
    let r = ((y_max - y) / raster.cell_size_y).floor();
    if c < 0.0 || r < 0.0 {
        return None;
    }
    let (c, r) = (c as usize, r as usize);
    (r < raster.rows && c < raster.cols).then_some((r, c))
}

/// Discards `pct` percent from each tail of a sorted slice.
fn trim_tails(sorted: &[f64], pct: f64) -> Vec<f64> {
    if pct <= 0.0 || sorted.len() < 3 {
        return sorted.to_vec();
    }
    let k = ((sorted.len() as f64) * pct / 100.0).floor() as usize;
    // Never trim away everything: with a small pool the tails would meet.
    let k = k.min((sorted.len() - 1) / 2);
    sorted[k..sorted.len() - k].to_vec()
}

fn reduce(sorted: &[f64], statistic: &str) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    match statistic {
        "median" => {
            let n = sorted.len();
            if n % 2 == 1 {
                sorted[n / 2]
            } else {
                (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
            }
        }
        "max" => sorted[sorted.len() - 1],
        "min" => sorted[0],
        _ => sorted.iter().sum::<f64>() / sorted.len() as f64,
    }
}

/// Sample standard deviation; zero for a single observation.
fn stddev(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (v.len() - 1) as f64;
    var.sqrt()
}

fn parse_wavelengths(
    args: &ToolArgs,
    key: &str,
    bands: usize,
) -> Result<Option<Vec<f64>>, ToolError> {
    let Some(s) = args.get(key).and_then(Value::as_str) else {
        return Ok(None);
    };
    if s.trim().is_empty() {
        return Ok(None);
    }
    let vals: Vec<f64> = s
        .split([',', ';'])
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(|t| {
            t.parse::<f64>()
                .map_err(|_| ToolError::Validation(format!("'{key}' entry '{t}' is not a number")))
        })
        .collect::<Result<_, _>>()?;
    if vals.len() != bands {
        return Err(ToolError::Validation(format!(
            "'{key}' has {} value(s) but the image has {bands} band(s)",
            vals.len()
        )));
    }
    Ok(Some(vals))
}

fn label_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.trim().to_string(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => fmt(f),
        FieldValue::Boolean(b) => b.to_string(),
        _ => String::new(),
    }
}

/// The library parser splits on commas and has no quoting, so a comma inside a
/// class name would shift every band value by one column. Substituting is the
/// only safe option and is cheaper than failing on a legitimate label.
fn csv_escape(s: &str) -> String {
    s.replace([',', '\n', '\r'], "_")
}

fn fmt(v: &f64) -> String {
    if v.is_finite() {
        let s = format!("{v:.6}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        "0".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};
    use wbvector::{Coord, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn tmp(tag: &str) -> String {
        static N: AtomicU64 = AtomicU64::new(0);
        let n = N.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!("esfi_{tag}_{}_{n}.csv", std::process::id()))
            .to_string_lossy()
            .to_string()
    }

    /// 4x4 image, `bands` bands. Band b at cell (r,c) = base(r,c) + 10*b.
    fn image(bands: usize) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 4,
            rows: 4,
            bands,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for b in 0..bands {
            for row in 0..4 {
                for col in 0..4 {
                    // Left half = 1.0, right half = 5.0, plus the band offset.
                    let base = if col < 2 { 1.0 } else { 5.0 };
                    r.set(
                        b as isize,
                        row as isize,
                        col as isize,
                        base + 10.0 * b as f64,
                    )
                    .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// Two polygons: "left" over columns 0-1, "right" over columns 2-3.
    fn training() -> String {
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        for (name, x0, x1) in [("left", 0.0, 2.0), ("right", 2.0, 4.0)] {
            l.add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(x0, 0.0),
                        Coord::xy(x1, 0.0),
                        Coord::xy(x1, 4.0),
                        Coord::xy(x0, 4.0),
                        Coord::xy(x0, 0.0),
                    ],
                    Vec::new(),
                )),
                &[("cls", name.into())],
            )
            .unwrap();
        }
        store(l)
    }

    fn store(l: Layer) -> String {
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (String, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ExtractSpectraFromImageTool.run(&args, &ctx()).unwrap();
        let p = res.outputs["output"].as_str().unwrap();
        let text = std::fs::read_to_string(p).unwrap();
        let _ = std::fs::remove_file(p);
        (text, res)
    }

    #[test]
    fn writes_one_line_per_class_with_one_value_per_band() {
        let (text, res) = run(json!({
            "input": image(3), "training_features": training(),
            "class_field": "cls", "output": tmp("basic"),
        }));
        let lines: Vec<&str> = text.lines().filter(|l| !l.is_empty()).collect();
        assert_eq!(lines.len(), 2);
        for l in &lines {
            assert_eq!(l.split(',').count(), 4, "name + 3 bands: {l}");
        }
        assert_eq!(res.outputs["class_count"], json!(2));
        assert_eq!(res.outputs["band_count"], json!(3));
    }

    #[test]
    fn the_spectra_are_the_class_means_per_band() {
        let (text, _) = run(json!({
            "input": image(3), "training_features": training(),
            "class_field": "cls", "output": tmp("means"),
        }));
        let left = text.lines().find(|l| l.starts_with("left,")).unwrap();
        // Left half is 1.0 in band 1, 11.0 in band 2, 21.0 in band 3.
        assert_eq!(left, "left,1,11,21", "got: {left}");
        let right = text.lines().find(|l| l.starts_with("right,")).unwrap();
        assert_eq!(right, "right,5,15,25", "got: {right}");
    }

    #[test]
    fn the_output_is_readable_by_spectral_library_matching() {
        // The requirement that makes the tool worth having: the exact CSV
        // layout the bundled consumer parses — headerless 'name,b1,b2,...'
        // with exactly one value per band and '#' comments allowed.
        let (text, _) = run(json!({
            "input": image(3), "training_features": training(),
            "class_field": "cls", "output": tmp("fmt"),
            "wavelengths": "450,550,650",
        }));
        let mut saw_comment = false;
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            if line.starts_with('#') {
                saw_comment = true;
                continue;
            }
            let parts: Vec<&str> = line.split(',').collect();
            assert!(
                parts[0].parse::<f64>().is_err(),
                "the first field must be a non-numeric name: {line}"
            );
            assert_eq!(parts.len() - 1, 3, "exactly one value per band: {line}");
            for p in &parts[1..] {
                assert!(p.parse::<f64>().is_ok(), "band value not numeric: {p}");
            }
        }
        assert!(saw_comment, "wavelengths should be recorded as a # comment");
    }

    #[test]
    fn statistics_other_than_mean_are_honoured() {
        // A polygon spanning both halves: mean 3, median 3, min 1, max 5.
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(4.0, 0.0),
                    Coord::xy(4.0, 4.0),
                    Coord::xy(0.0, 4.0),
                    Coord::xy(0.0, 0.0),
                ],
                Vec::new(),
            )),
            &[("cls", "all".into())],
        )
        .unwrap();
        let path = store(l);
        for (stat, expect) in [("mean", "3"), ("median", "3"), ("min", "1"), ("max", "5")] {
            let (text, _) = run(json!({
                "input": image(1), "training_features": path,
                "class_field": "cls", "output": tmp(stat), "statistic": stat,
            }));
            assert_eq!(text.trim(), format!("all,{expect}"), "statistic {stat}");
        }
    }

    #[test]
    fn a_cell_with_any_no_data_band_is_excluded_entirely() {
        // Otherwise the bands of one "spectrum" would come from different cells.
        let mut r = Raster::new(RasterConfig {
            cols: 2,
            rows: 1,
            bands: 2,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        // Cell 0 is complete; cell 1 is no-data in band 2 only.
        r.set(0, 0, 0, 2.0).unwrap();
        r.set(1, 0, 0, 4.0).unwrap();
        r.set(0, 0, 1, 100.0).unwrap();
        r.set(1, 0, 1, -9999.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let img = wbraster::memory_store::make_raster_memory_path(&id);

        let mut l = Layer::new("t").with_geom_type(GeometryType::Polygon);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(2.0, 0.0),
                    Coord::xy(2.0, 1.0),
                    Coord::xy(0.0, 1.0),
                    Coord::xy(0.0, 0.0),
                ],
                Vec::new(),
            )),
            &[("cls", "c".into())],
        )
        .unwrap();
        let (text, _) = run(json!({
            "input": img, "training_features": store(l),
            "class_field": "cls", "output": tmp("nodata"),
        }));
        // Only cell 0 contributes: the 100.0 in band 1 must not leak in.
        assert_eq!(text.trim(), "c,2,4");
    }

    #[test]
    fn trimming_suppresses_a_bright_outlier() {
        // One bright roof pixel inside a vegetation polygon is exactly what
        // shifts an endmember enough to mis-classify a scene.
        let mut r = Raster::new(RasterConfig {
            cols: 10,
            rows: 1,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for c in 0..9 {
            r.set(0, 0, c, 10.0).unwrap();
        }
        r.set(0, 0, 9, 1000.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let img = wbraster::memory_store::make_raster_memory_path(&id);

        let mut l = Layer::new("t").with_geom_type(GeometryType::Polygon);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(10.0, 0.0),
                    Coord::xy(10.0, 1.0),
                    Coord::xy(0.0, 1.0),
                    Coord::xy(0.0, 0.0),
                ],
                Vec::new(),
            )),
            &[("cls", "veg".into())],
        )
        .unwrap();
        let path = store(l);

        let (untrimmed, _) = run(json!({
            "input": img, "training_features": path,
            "class_field": "cls", "output": tmp("t0"),
        }));
        assert_eq!(untrimmed.trim(), "veg,109", "outlier dominates the mean");

        let (trimmed, _) = run(json!({
            "input": img, "training_features": path,
            "class_field": "cls", "output": tmp("t1"), "trim_percent": 10.0,
        }));
        assert_eq!(trimmed.trim(), "veg,10", "trimming removes it");
    }

    #[test]
    fn point_training_features_sample_their_own_cell() {
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_feature(Some(Geometry::point(0.5, 0.5)), &[("cls", "a".into())])
            .unwrap();
        l.add_feature(Some(Geometry::point(3.5, 0.5)), &[("cls", "b".into())])
            .unwrap();
        let (text, res) = run(json!({
            "input": image(1), "training_features": store(l),
            "class_field": "cls", "output": tmp("pts"),
        }));
        assert_eq!(res.outputs["class_count"], json!(2));
        assert!(text.contains("a,1"), "got: {text}");
        assert!(text.contains("b,5"), "got: {text}");
    }

    #[test]
    fn features_sharing_a_label_pool_into_one_spectrum() {
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        for x in [0.5, 1.5, 2.5, 3.5] {
            l.add_feature(Some(Geometry::point(x, 0.5)), &[("cls", "all".into())])
                .unwrap();
        }
        let (text, res) = run(json!({
            "input": image(1), "training_features": store(l),
            "class_field": "cls", "output": tmp("pool"),
        }));
        assert_eq!(res.outputs["class_count"], json!(1));
        // Mean of 1,1,5,5.
        assert_eq!(text.trim(), "all,3");
    }

    #[test]
    fn a_class_below_min_cells_is_skipped_and_reported() {
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_feature(Some(Geometry::point(0.5, 0.5)), &[("cls", "rare".into())])
            .unwrap();
        for x in [2.5, 3.5] {
            l.add_feature(Some(Geometry::point(x, 0.5)), &[("cls", "big".into())])
                .unwrap();
        }
        let (text, res) = run(json!({
            "input": image(1), "training_features": store(l),
            "class_field": "cls", "output": tmp("min"), "min_cells": 2,
        }));
        assert_eq!(res.outputs["skipped_classes"], json!(["rare"]));
        assert!(!text.contains("rare"));
        assert!(text.contains("big"));
    }

    #[test]
    fn the_stats_file_carries_the_spread_and_counts() {
        let stats = tmp("stats");
        let (_, res) = run(json!({
            "input": image(2), "training_features": training(),
            "class_field": "cls", "output": tmp("lib"),
            "output_stats": stats.clone(), "wavelengths": "450,550",
        }));
        assert_eq!(res.outputs["output_stats"], json!(stats));
        let text = std::fs::read_to_string(&stats).unwrap();
        let _ = std::fs::remove_file(&stats);
        assert!(text.starts_with("class,band,wavelength,value,stddev,count"));
        // Each class is uniform within its half, so the spread is zero and
        // each class covers 8 of the 16 cells.
        assert!(text.contains("left,1,450,1,0,8"), "got: {text}");
    }

    #[test]
    fn a_comma_in_a_class_name_cannot_shift_the_band_columns() {
        // The library parser has no quoting, so an unescaped comma would push
        // every band value one column to the right.
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_feature(
            Some(Geometry::point(0.5, 0.5)),
            &[("cls", "urban, dense".into())],
        )
        .unwrap();
        let (text, _) = run(json!({
            "input": image(2), "training_features": store(l),
            "class_field": "cls", "output": tmp("comma"),
        }));
        let line = text.lines().find(|l| !l.starts_with('#')).unwrap();
        assert_eq!(line.split(',').count(), 3, "name + 2 bands: {line}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExtractSpectraFromImageTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "i.tif"})));
        let base = json!({
            "input": "i.tif", "training_features": "t.shp",
            "class_field": "c", "output": "o.csv",
        });
        let with = |k: &str, v: Value| {
            let mut m = base.clone();
            m[k] = v;
            m
        };
        assert!(bad(with("statistic", json!("mode"))));
        // 50% from each tail would discard everything.
        assert!(bad(with("trim_percent", json!(50))));
        assert!(bad(with("trim_percent", json!(-1))));
        assert!(bad(with("min_cells", json!(0))));
    }
}
