//! Shared loading and validation for co-registered raster *stacks*.
//!
//! `cell_statistics` established the stack contract for GeoLibre's local
//! (multi-raster, per-cell) tools: inputs are one multiband raster or a
//! comma-separated list of rasters, every band of every raster becomes a layer,
//! and the rasters must share a grid. `cell_position_statistics` and
//! `frequency_comparison` need exactly the same contract, so it lives here
//! rather than being copied a third time.
//!
//! Two band policies are supported, matching ArcGIS's `SINGLE_BAND` /
//! `MULTI_BAND` process-as flag:
//!
//! * [`BandPolicy::SingleBand`] — every band of every input is one layer of a
//!   single stack, producing a single-band result. This is what
//!   `cell_statistics` does.
//! * [`BandPolicy::MultiBand`] — band *i* of every input forms its own stack,
//!   producing an *n*-band result. Requires all inputs to have the same band
//!   count, which is checked up front rather than discovered mid-run.

use serde_json::Value;
use wbcore::{ToolArgs, ToolError};
use wbraster::{DataType, Raster, RasterConfig};

/// How input bands are grouped into stacks.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum BandPolicy {
    /// All bands of all inputs form one stack → one output band.
    SingleBand,
    /// Band *i* of every input forms stack *i* → *n* output bands.
    MultiBand,
}

impl BandPolicy {
    pub(crate) fn label(self) -> &'static str {
        match self {
            BandPolicy::SingleBand => "single_band",
            BandPolicy::MultiBand => "multi_band",
        }
    }
}

/// A co-registered set of rasters, pre-grouped into per-output-band stacks.
pub(crate) struct Stack {
    rasters: Vec<Raster>,
    /// One entry per output band; each is the list of `(raster index, band)`
    /// layers that reduce into that output band.
    groups: Vec<Vec<(usize, isize)>>,
    pub(crate) rows: usize,
    pub(crate) cols: usize,
}

impl Stack {
    /// Number of output bands.
    pub(crate) fn group_count(&self) -> usize {
        self.groups.len()
    }

    /// Total input layers across all groups (for reporting).
    pub(crate) fn total_layers(&self) -> usize {
        self.groups.iter().map(Vec::len).sum()
    }

    /// The first raster, used as the output geometry/CRS template.
    pub(crate) fn template(&self) -> &Raster {
        &self.rasters[0]
    }

    /// Collects the valid observations at `(row, col)` for output band `group`
    /// into `out`, returning `true` if any layer was no-data.
    ///
    /// `out` is cleared first, so callers can reuse one buffer across cells.
    pub(crate) fn cell_values(
        &self,
        group: usize,
        row: usize,
        col: usize,
        out: &mut Vec<f64>,
    ) -> bool {
        out.clear();
        let mut had_nodata = false;
        for &(ri, band) in &self.groups[group] {
            let ras = &self.rasters[ri];
            let v = ras.get(band, row as isize, col as isize);
            if v != ras.nodata && v.is_finite() {
                out.push(v);
            } else {
                had_nodata = true;
            }
        }
        had_nodata
    }
}

/// Splits a comma-separated `inputs` parameter into paths.
pub(crate) fn parse_input_paths(args: &ToolArgs, key: &str) -> Result<Vec<String>, ToolError> {
    let s = args
        .get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))?;
    let paths: Vec<String> = s
        .split(',')
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(String::from)
        .collect();
    if paths.is_empty() {
        return Err(ToolError::Validation(format!("'{key}' is empty")));
    }
    Ok(paths)
}

/// Parses the ArcGIS-style process-as flag.
pub(crate) fn parse_band_policy(args: &ToolArgs, key: &str) -> Result<BandPolicy, ToolError> {
    match args.get(key).and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("single_band") | Some("SINGLE_BAND") => Ok(BandPolicy::SingleBand),
        Some("multi_band") | Some("MULTI_BAND") => Ok(BandPolicy::MultiBand),
        Some(o) => Err(ToolError::Validation(format!(
            "'{key}' must be 'single_band' or 'multi_band', got '{o}'"
        ))),
    }
}

/// Loads and co-registration-checks a raster stack.
///
/// `min_layers` is the smallest number of layers each group must contain; a
/// stack reducer needs at least 2, while a comparison against a separate
/// reference raster is meaningful with 1.
pub(crate) fn load_stack(
    paths: &[String],
    policy: BandPolicy,
    min_layers: usize,
) -> Result<Stack, ToolError> {
    if paths.is_empty() {
        return Err(ToolError::Validation(
            "no input rasters were supplied".to_string(),
        ));
    }
    let rasters: Vec<Raster> = paths
        .iter()
        .map(|p| crate::common::load_input_raster(p))
        .collect::<Result<_, _>>()?;
    let base_rows = rasters[0].rows;
    let base_cols = rasters[0].cols;
    check_alignment(&rasters)?;

    let groups: Vec<Vec<(usize, isize)>> = match policy {
        BandPolicy::SingleBand => {
            let mut layers = Vec::new();
            for (i, r) in rasters.iter().enumerate() {
                for b in 0..r.bands {
                    layers.push((i, b as isize));
                }
            }
            vec![layers]
        }
        BandPolicy::MultiBand => {
            let bands = rasters[0].bands;
            for (i, r) in rasters.iter().enumerate() {
                if r.bands != bands {
                    return Err(ToolError::Validation(format!(
                        "process_as_multiband requires every input to have the same band count; \
                         input 0 has {bands}, input {i} has {}",
                        r.bands
                    )));
                }
            }
            (0..bands)
                .map(|b| (0..rasters.len()).map(|i| (i, b as isize)).collect())
                .collect()
        }
    };

    for (g, layers) in groups.iter().enumerate() {
        if layers.len() < min_layers {
            return Err(ToolError::Validation(format!(
                "band group {g} has {} layer(s); need at least {min_layers}",
                layers.len()
            )));
        }
    }

    Ok(Stack {
        rasters,
        groups,
        rows: base_rows,
        cols: base_cols,
    })
}

/// Rejects rasters that are not on the same grid.
///
/// Cell values are combined position-by-position, so same-size rasters from
/// different geotransforms would otherwise be silently blended under the first
/// raster's geometry — the check `cell_statistics` already performs.
///
/// Takes references: the check reads only headers, and an owning signature
/// would force every ad-hoc two-raster caller to `clone()` both operands,
/// doubling peak memory exactly where the tools allocate their working buffers.
pub(crate) fn check_alignment_refs(rasters: &[&Raster]) -> Result<(), ToolError> {
    let Some(base) = rasters.first() else {
        return Ok(());
    };
    let (rows, cols) = (base.rows, base.cols);
    let aligned = |a: f64, b: f64| (a - b).abs() <= 1e-6 * a.abs().max(b.abs()).max(1.0);
    for (i, r) in rasters.iter().enumerate().skip(1) {
        if r.rows != rows || r.cols != cols {
            return Err(ToolError::Validation(format!(
                "raster {i} is {}x{}, expected {rows}x{cols}",
                r.rows, r.cols
            )));
        }
        if !aligned(r.x_min, base.x_min)
            || !aligned(r.y_min, base.y_min)
            || !aligned(r.cell_size_x, base.cell_size_x)
            || !aligned(r.cell_size_y, base.cell_size_y)
        {
            return Err(ToolError::Validation(format!(
                "raster {i} is not co-registered with input 0 (origin/resolution differ); \
                 inputs must share the same grid"
            )));
        }
        if r.crs.epsg.is_some() && base.crs.epsg.is_some() && r.crs.epsg != base.crs.epsg {
            return Err(ToolError::Validation(format!(
                "raster {i} CRS (EPSG {:?}) differs from input 0 (EPSG {:?})",
                r.crs.epsg, base.crs.epsg
            )));
        }
    }
    Ok(())
}

/// Owning-slice wrapper over [`check_alignment_refs`], for callers that already
/// hold a `Vec<Raster>`.
pub(crate) fn check_alignment(rasters: &[Raster]) -> Result<(), ToolError> {
    check_alignment_refs(&rasters.iter().collect::<Vec<_>>())
}

/// Builds an output raster with `bands` bands from `template`'s geometry.
///
/// `common::raster_like_with_data` allocates a single band, which is the right
/// default for most tools; the `MULTI_BAND` policy needs an *n*-band store, so
/// this is its multiband sibling. `data` is band-major: band *b* occupies
/// `data[b][row * cols + col]`.
pub(crate) fn raster_like_multiband(
    template: &Raster,
    data: &[Vec<f64>],
    nodata: f64,
    data_type: DataType,
) -> Result<Raster, ToolError> {
    let rows = template.rows;
    let cols = template.cols;
    let mut out = Raster::new(RasterConfig {
        cols,
        rows,
        bands: data.len(),
        x_min: template.x_min,
        y_min: template.y_min,
        cell_size: template.cell_size_x,
        cell_size_y: Some(template.cell_size_y),
        nodata,
        data_type,
        crs: template.crs.clone(),
        metadata: template.metadata.clone(),
    });
    for (b, band) in data.iter().enumerate() {
        if band.len() != rows * cols {
            return Err(ToolError::Execution(format!(
                "output band {b} length {} does not match {rows}x{cols}",
                band.len()
            )));
        }
        for row in 0..rows {
            for col in 0..cols {
                out.set(
                    b as isize,
                    row as isize,
                    col as isize,
                    band[row * cols + col],
                )
                .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
            }
        }
    }
    Ok(out)
}

/// Writes a stack result, choosing the single- or multi-band writer.
pub(crate) fn write_stack_result(
    template: &Raster,
    bands: Vec<Vec<f64>>,
    nodata: f64,
    data_type: DataType,
    output: Option<&str>,
) -> Result<String, ToolError> {
    let raster = if bands.len() == 1 {
        crate::common::raster_like_with_data(
            template,
            bands.into_iter().next().expect("one band"),
            nodata,
            data_type,
        )?
    } else {
        raster_like_multiband(template, &bands, nodata, data_type)?
    };
    crate::common::write_or_store_output(raster, output)
}
