//! Converts the GeoJSON geometry objects `polygonize_to_geojson` emits back
//! into `wbvector::Geometry`.
//!
//! `polygonize` returns a GeoJSON string because that is what the GeoLibre /
//! MapLibre side consumes directly, so every tool that wants the traced rings
//! as features has to parse them again. `percentile_contours` grew a private
//! copy of this first; the round-18 mask tools need the same thing, so it lives
//! here.

use serde_json::Value;
use wbvector::{Coord, Geometry, Ring};

/// Parses one GeoJSON `Polygon` object (`{"type": "Polygon", "coordinates":
/// [...]}`) into a [`Geometry::Polygon`].
///
/// Returns `None` when the object is not a polygon or its exterior ring has
/// fewer than three distinct vertices. Rings arrive closed (first point
/// repeated) and are stored unclosed, which is the crate's convention.
pub(crate) fn geometry_from_json(geom: &Value) -> Option<Geometry> {
    if geom.get("type").and_then(Value::as_str) != Some("Polygon") {
        return None;
    }
    let coords = geom.get("coordinates")?.as_array()?;
    let mut rings = coords.iter().filter_map(ring_from_json);
    let exterior = rings.next()?;
    let interiors: Vec<Ring> = rings.collect();
    Some(Geometry::Polygon {
        exterior,
        interiors,
    })
}

/// Parses one GeoJSON linear ring, dropping the repeated closing vertex.
pub(crate) fn ring_from_json(ring: &Value) -> Option<Ring> {
    let pts = ring.as_array()?;
    let mut coords: Vec<Coord> = pts
        .iter()
        .filter_map(|p| {
            let a = p.as_array()?;
            Some(Coord::xy(a.first()?.as_f64()?, a.get(1)?.as_f64()?))
        })
        .collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    if coords.len() < 3 {
        return None;
    }
    Some(Ring::new(coords))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn parses_polygon_with_hole() {
        let g = json!({
            "type": "Polygon",
            "coordinates": [
                [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0], [0.0, 0.0]],
                [[1.0, 1.0], [1.0, 2.0], [2.0, 2.0], [2.0, 1.0], [1.0, 1.0]]
            ]
        });
        let Some(Geometry::Polygon {
            exterior,
            interiors,
        }) = geometry_from_json(&g)
        else {
            panic!("expected a polygon");
        };
        // Closing vertex dropped on both rings.
        assert_eq!(exterior.coords().len(), 4);
        assert_eq!(interiors.len(), 1);
        assert_eq!(interiors[0].coords().len(), 4);
    }

    #[test]
    fn rejects_non_polygons_and_degenerate_rings() {
        assert!(geometry_from_json(&json!({"type": "Point", "coordinates": [0, 0]})).is_none());
        let degenerate = json!({
            "type": "Polygon",
            "coordinates": [[[0.0, 0.0], [1.0, 1.0], [0.0, 0.0]]]
        });
        assert!(geometry_from_json(&degenerate).is_none());
    }
}
