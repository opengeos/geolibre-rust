//! Shared raster rendering: map raster values through a colormap to RGBA and
//! encode PNGs. Used by `render_raster_png` (one image) and `raster_to_tiles`
//! (a pyramid of 256x256 web-map tiles).
//!
//! Colormaps are stored as a handful of evenly spaced RGB control stops and
//! linearly interpolated, which keeps the table tiny while looking smooth.

use wbcore::ToolError;

/// A perceptual or terrain colormap selectable by name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Colormap {
    Grayscale,
    Viridis,
    Magma,
    Terrain,
    Turbo,
}

impl Colormap {
    /// Parses a colormap name (case-insensitive). Defaults are chosen by the
    /// caller; this returns an error for unknown names so typos are caught.
    pub fn parse(name: &str) -> Result<Self, ToolError> {
        match name.trim().to_ascii_lowercase().as_str() {
            "grayscale" | "greyscale" | "gray" | "grey" => Ok(Self::Grayscale),
            "viridis" => Ok(Self::Viridis),
            "magma" => Ok(Self::Magma),
            "terrain" => Ok(Self::Terrain),
            "turbo" => Ok(Self::Turbo),
            other => Err(ToolError::Validation(format!(
                "unknown colormap '{other}' (expected grayscale, viridis, magma, terrain, or turbo)"
            ))),
        }
    }

    /// Evenly spaced RGB control stops for this colormap.
    fn stops(self) -> &'static [[u8; 3]] {
        match self {
            Self::Grayscale => &[[0, 0, 0], [255, 255, 255]],
            // Compact subsamples of matplotlib's viridis / magma / turbo.
            Self::Viridis => &[
                [68, 1, 84],
                [72, 40, 120],
                [62, 74, 137],
                [49, 104, 142],
                [38, 130, 142],
                [31, 158, 137],
                [53, 183, 121],
                [110, 206, 88],
                [181, 222, 43],
                [253, 231, 37],
            ],
            Self::Magma => &[
                [0, 0, 4],
                [28, 16, 68],
                [79, 18, 123],
                [129, 37, 129],
                [181, 54, 122],
                [229, 80, 100],
                [251, 135, 97],
                [254, 194, 135],
                [252, 253, 191],
            ],
            Self::Turbo => &[
                [48, 18, 59],
                [70, 107, 227],
                [40, 187, 226],
                [49, 242, 153],
                [150, 254, 70],
                [225, 220, 55],
                [253, 152, 39],
                [223, 67, 19],
                [122, 4, 3],
            ],
            // A hypsometric land ramp: greens -> tans -> browns -> white.
            Self::Terrain => &[
                [0, 97, 71],
                [54, 145, 79],
                [149, 191, 99],
                [223, 217, 140],
                [199, 165, 110],
                [151, 110, 75],
                [132, 95, 90],
                [255, 255, 255],
            ],
        }
    }

    /// Maps a normalized value in [0, 1] to an RGB triple via piecewise-linear
    /// interpolation between control stops. Values are clamped to [0, 1].
    pub fn rgb(self, t: f64) -> [u8; 3] {
        let stops = self.stops();
        let t = t.clamp(0.0, 1.0);
        if stops.len() == 1 {
            return stops[0];
        }
        let scaled = t * (stops.len() - 1) as f64;
        let i = (scaled.floor() as usize).min(stops.len() - 2);
        let frac = scaled - i as f64;
        let a = stops[i];
        let b = stops[i + 1];
        [
            lerp(a[0], b[0], frac),
            lerp(a[1], b[1], frac),
            lerp(a[2], b[2], frac),
        ]
    }
}

fn lerp(a: u8, b: u8, frac: f64) -> u8 {
    (a as f64 + (b as f64 - a as f64) * frac).round().clamp(0.0, 255.0) as u8
}

/// Linear stretch from a value range to [0, 1]. A zero-width range maps every
/// valid value to 0.0 so a constant raster still renders.
#[inline]
pub fn normalize(value: f64, min: f64, max: f64) -> f64 {
    if max <= min {
        0.0
    } else {
        (value - min) / (max - min)
    }
}

/// Encodes an RGBA8 buffer (`width*height*4` bytes, row-major) as a PNG byte
/// stream. Alpha 0 marks no-data / out-of-bounds pixels, so tiles overlay
/// cleanly on a web map.
pub fn encode_png_rgba(rgba: &[u8], width: u32, height: u32) -> Result<Vec<u8>, ToolError> {
    let expected = width as usize * height as usize * 4;
    if rgba.len() != expected {
        return Err(ToolError::Execution(format!(
            "RGBA buffer length {} does not match {width}x{height}x4 = {expected}",
            rgba.len()
        )));
    }
    let mut out = Vec::new();
    {
        let mut encoder = png::Encoder::new(&mut out, width, height);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        let mut writer = encoder
            .write_header()
            .map_err(|e| ToolError::Execution(format!("failed writing PNG header: {e}")))?;
        writer
            .write_image_data(rgba)
            .map_err(|e| ToolError::Execution(format!("failed writing PNG data: {e}")))?;
    }
    Ok(out)
}
