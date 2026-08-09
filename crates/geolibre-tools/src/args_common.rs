//! Small shared argument parsers for GeoLibre tools.
//!
//! `common.rs` covers raster I/O and `vector_common.rs` covers vector I/O, but
//! the scalar parsers (`bool`, `f64`, `usize`, required strings) have
//! historically been redeclared privately in each tool module. That is fine for
//! one tool and wasteful across a batch, so the canonical versions live here.
//!
//! Semantics match the existing private copies exactly, so behaviour is
//! unchanged for callers: JSON booleans/numbers pass through, strings are
//! parsed, and an empty string means "not supplied" (the same convention
//! `common::parse_optional_output` uses when it maps `""` to `None`).

use serde_json::Value;
use wbcore::{ToolArgs, ToolError};

/// Optional boolean. Accepts JSON `true`/`false`, numbers (non-zero = true),
/// and the strings `true`/`1`/`yes` and `false`/`0`/`no` (case-insensitive).
pub(crate) fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(0.0) != 0.0)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

/// Optional boolean with a default.
pub(crate) fn bool_or(args: &ToolArgs, key: &str, default: bool) -> Result<bool, ToolError> {
    Ok(opt_bool(args, key)?.unwrap_or(default))
}

/// Optional finite `f64`. Rejects NaN/infinity so downstream geometry never
/// silently propagates a non-finite coordinate.
pub(crate) fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    let v = match args.get(key) {
        None | Some(Value::Null) => return Ok(None),
        Some(Value::Number(n)) => n.as_f64().unwrap_or(f64::NAN),
        Some(Value::String(s)) if s.trim().is_empty() => return Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number")))?,
        Some(_) => {
            return Err(ToolError::Validation(format!(
                "parameter '{key}' must be a number"
            )))
        }
    };
    if !v.is_finite() {
        return Err(ToolError::Validation(format!(
            "parameter '{key}' must be a finite number"
        )));
    }
    Ok(Some(v))
}

/// Optional `f64` with a default.
pub(crate) fn f64_or(args: &ToolArgs, key: &str, default: f64) -> Result<f64, ToolError> {
    Ok(opt_f64(args, key)?.unwrap_or(default))
}

/// Optional strictly-positive `f64`.
pub(crate) fn opt_positive_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        None => Ok(None),
        Some(v) if v > 0.0 => Ok(Some(v)),
        Some(v) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be > 0, got {v}"
        ))),
    }
}

/// Optional non-negative integer.
pub(crate) fn opt_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match opt_f64(args, key)? {
        None => Ok(None),
        Some(v) if v >= 0.0 && v.fract() == 0.0 => Ok(Some(v as usize)),
        Some(v) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a non-negative whole number, got {v}"
        ))),
    }
}

/// Optional integer with a default.
pub(crate) fn usize_or(args: &ToolArgs, key: &str, default: usize) -> Result<usize, ToolError> {
    Ok(opt_usize(args, key)?.unwrap_or(default))
}

/// Required non-empty string.
pub(crate) fn req_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    let s = args
        .get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))?;
    Ok(s)
}

/// Optional lower-cased enum-ish string, with `""` treated as absent.
pub(crate) fn opt_choice(args: &ToolArgs, key: &str) -> Option<String> {
    args.get(key)
        .and_then(Value::as_str)
        .map(|s| s.trim().to_ascii_lowercase())
        .filter(|s| !s.is_empty())
}

/// Resolves a choice parameter against an allowed set, defaulting when absent.
pub(crate) fn choice_or<'a>(
    args: &ToolArgs,
    key: &str,
    allowed: &[&'a str],
    default: &'a str,
) -> Result<&'a str, ToolError> {
    match opt_choice(args, key) {
        None => Ok(default),
        Some(v) => allowed
            .iter()
            .find(|a| **a == v)
            .copied()
            .ok_or_else(|| {
                ToolError::Validation(format!(
                    "'{key}' must be one of {}, got '{v}'",
                    allowed.join("|")
                ))
            }),
    }
}

/// 1-based band selector, returned as the 0-based index `wbraster` expects.
pub(crate) fn band_index(args: &ToolArgs, key: &str) -> Result<isize, ToolError> {
    match opt_usize(args, key)? {
        None => Ok(0),
        Some(0) => Err(ToolError::Validation(format!(
            "'{key}' is 1-based; use 1 for the first band"
        ))),
        Some(b) => Ok(b as isize - 1),
    }
}
