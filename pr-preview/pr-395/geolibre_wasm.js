/* @ts-self-types="./geolibre_wasm.d.ts" */

/**
 * Builder for encoding a Cloud Optimized GeoTIFF (tiled, with overviews and
 * GDAL ghost metadata) to bytes. A COG is also a valid plain GeoTIFF.
 *
 * Configure with the `set_*` methods, then call one of `write_*` with the
 * pixel data to get a `Uint8Array` of the encoded file.
 */
export class CogBuilder {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CogBuilderFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_cogbuilder_free(ptr, 0);
    }
    /**
     * New builder for a `width` x `height` raster with `bands` bands.
     * @param {number} width
     * @param {number} height
     * @param {number} bands
     */
    constructor(width, height, bands) {
        const ret = wasm.cogbuilder_new(width, height, bands);
        this.__wbg_ptr = ret;
        CogBuilderFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Force BigTIFF (64-bit offsets) for very large outputs.
     * @param {boolean} on
     */
    set_bigtiff(on) {
        wasm.cogbuilder_set_bigtiff(this.__wbg_ptr, on);
    }
    /**
     * Compression: `none`, `lzw`, `deflate`, `packbits`, `webp`, `jpeg`, `jpegxl`.
     * @param {string} name
     */
    set_compression(name) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passStringToWasm0(name, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogbuilder_set_compression(retptr, this.__wbg_ptr, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            if (r1) {
                throw takeObject(r0);
            }
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Set the EPSG code (1..=65535).
     * @param {number} epsg
     */
    set_epsg(epsg) {
        wasm.cogbuilder_set_epsg(this.__wbg_ptr, epsg);
    }
    /**
     * Set the full affine geo-transform:
     * `[x_origin, pixel_width, row_rotation, y_origin, col_rotation, pixel_height]`.
     * @param {Float64Array} gt
     */
    set_geo_transform(gt) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArrayF64ToWasm0(gt, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogbuilder_set_geo_transform(retptr, this.__wbg_ptr, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            if (r1) {
                throw takeObject(r0);
            }
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Set the no-data sentinel value.
     * @param {number} v
     */
    set_nodata(v) {
        wasm.cogbuilder_set_nodata(this.__wbg_ptr, v);
    }
    /**
     * Convenience: north-up geo-transform from upper-left origin and pixel size.
     * @param {number} x_min
     * @param {number} y_max
     * @param {number} pixel_size
     */
    set_origin(x_min, y_max, pixel_size) {
        wasm.cogbuilder_set_origin(this.__wbg_ptr, x_min, y_max, pixel_size);
    }
    /**
     * Explicit overview decimation factors (e.g. `[2,4,8]`); empty disables overviews.
     * @param {Uint32Array} levels
     */
    set_overview_levels(levels) {
        const ptr0 = passArray32ToWasm0(levels, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        wasm.cogbuilder_set_overview_levels(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * Internal tile size in pixels (default 512).
     * @param {number} px
     */
    set_tile_size(px) {
        wasm.cogbuilder_set_tile_size(this.__wbg_ptr, px);
    }
    /**
     * Encode `f32` pixel data to a COG. `Uint8Array`.
     * @param {Float32Array} data
     * @returns {Uint8Array}
     */
    write_f32(data) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArrayF32ToWasm0(data, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogbuilder_write_f32(retptr, this.__wbg_ptr, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v2 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 1, 1);
            return v2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Encode `f64` pixel data to a COG. `Uint8Array`.
     * @param {Float64Array} data
     * @returns {Uint8Array}
     */
    write_f64(data) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArrayF64ToWasm0(data, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogbuilder_write_f64(retptr, this.__wbg_ptr, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v2 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 1, 1);
            return v2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Encode `u8` pixel data to a COG. `Uint8Array`.
     * @param {Uint8Array} data
     * @returns {Uint8Array}
     */
    write_u8(data) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogbuilder_write_u8(retptr, this.__wbg_ptr, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v2 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 1, 1);
            return v2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
}
if (Symbol.dispose) CogBuilder.prototype[Symbol.dispose] = CogBuilder.prototype.free;

/**
 * Range-request reader for a (tiled) Cloud Optimized GeoTIFF.
 *
 * The wasm module does no network I/O itself; this class parses the header and
 * tells the JS host exactly which byte ranges to fetch, then decodes the tiles
 * the host fetches. Typical flow:
 *
 * 1. Range-fetch the first chunk of the file (e.g. 0..1 MiB) and
 *    `new CogStream(headerBytes)`. If it throws "need more header bytes", fetch
 *    a larger prefix and retry.
 * 2. Pick a level (0 = full res, higher = overviews) and a pixel window.
 * 3. `tiles_for_window(level, x, y, w, h)` returns the tiles and their byte
 *    ranges; range-fetch each, then `decode_tile_f64(level, bytes)`.
 */
export class CogStream {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CogStreamFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_cogstream_free(ptr, 0);
    }
    /**
     * Reproject a bbox from `bbox_epsg` into this COG's dataset CRS.
     *
     * The COG projection string is preferred over its EPSG tag when available,
     * because some user-defined projected GeoTIFFs expose only their geographic
     * base EPSG code.
     * @param {number} bbox_epsg
     * @param {Float64Array} bbox
     * @returns {Float64Array}
     */
    bbox_to_dataset_crs(bbox_epsg, bbox) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArrayF64ToWasm0(bbox, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogstream_bbox_to_dataset_crs(retptr, this.__wbg_ptr, bbox_epsg, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v2 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Bounding box `[min_x, min_y, max_x, max_y]` in the dataset CRS, or empty.
     * @returns {Float64Array}
     */
    bounding_box() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_bounding_box(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Bounds `[min_lon, min_lat, max_lon, max_lat]` in WGS84 degrees, or empty.
     * @returns {Float64Array}
     */
    bounds_lonlat() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_bounds_lonlat(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Image center `[x, y]` in the dataset CRS, or empty.
     * @returns {Float64Array}
     */
    center() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_center(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Image center `[lon, lat]` in WGS84 degrees, or empty if not convertible.
     * @returns {Float64Array}
     */
    center_lonlat() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_center_lonlat(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Decode one tile's fetched (compressed) bytes into an `f64` `Float64Array`,
     * pixel-interleaved, length `tile_width * tile_height * bands`. Edge tiles
     * come back full-size; clip to the image/window on the JS side.
     * @param {number} level
     * @param {Uint8Array} tile_bytes
     * @returns {Float64Array}
     */
    decode_tile_f64(level, tile_bytes) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArray8ToWasm0(tile_bytes, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogstream_decode_tile_f64(retptr, this.__wbg_ptr, level, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v2 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * EPSG code of the full-resolution level, if any.
     * @returns {number | undefined}
     */
    get epsg() {
        const ret = wasm.cogstream_epsg(this.__wbg_ptr);
        return ret === Number.MAX_SAFE_INTEGER ? undefined : ret;
    }
    /**
     * Level-0 geo-transform `[x_origin, pixel_width, row_rot, y_origin, col_rot,
     * pixel_height]`, or empty if not georeferenced.
     * @returns {Float64Array}
     */
    geo_transform() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_geo_transform(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * True when the COG CRS is represented by a user-defined projection string.
     * @returns {boolean}
     */
    get has_projection_string() {
        const ret = wasm.cogstream_has_projection_string(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * JSON array describing every level: `[{level,width,height,tile_width,
     * tile_height,tiles_x,tiles_y,bands,bits_per_sample,sample_format,compression}]`.
     * @returns {string}
     */
    levels_json() {
        let deferred1_0;
        let deferred1_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_levels_json(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            deferred1_0 = r0;
            deferred1_1 = r1;
            return getStringFromWasm0(r0, r1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Parse a COG's tile layout from front-of-file header bytes.
     * @param {Uint8Array} header_bytes
     */
    constructor(header_bytes) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArray8ToWasm0(header_bytes, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogstream_new(retptr, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            if (r2) {
                throw takeObject(r1);
            }
            this.__wbg_ptr = r0;
            CogStreamFinalization.register(this, this.__wbg_ptr, this);
            return this;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * No-data sentinel, if declared.
     * @returns {number | undefined}
     */
    get nodata() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_nodata(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r2 = getDataViewMemory0().getFloat64(retptr + 8 * 1, true);
            return r0 === 0 ? undefined : r2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Number of resolution levels (1 + overview count).
     * @returns {number}
     */
    get num_levels() {
        const ret = wasm.cogstream_num_levels(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Reproject x,y coordinate pairs from this COG's dataset CRS to an EPSG CRS.
     * @param {number} dst_epsg
     * @param {Float64Array} xy
     * @returns {Float64Array}
     */
    points_from_dataset_crs(dst_epsg, xy) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArrayF64ToWasm0(xy, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogstream_points_from_dataset_crs(retptr, this.__wbg_ptr, dst_epsg, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v2 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Reproject x,y coordinate pairs from an EPSG CRS into this COG's dataset CRS.
     * @param {number} src_epsg
     * @param {Float64Array} xy
     * @returns {Float64Array}
     */
    points_to_dataset_crs(src_epsg, xy) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArrayF64ToWasm0(xy, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.cogstream_points_to_dataset_crs(retptr, this.__wbg_ptr, src_epsg, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v2 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * `[offset, length]` byte range of the tile at `(col, row)` on `level`.
     * @param {number} level
     * @param {number} col
     * @param {number} row
     * @returns {Float64Array}
     */
    tile_range(level, col, row) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_tile_range(retptr, this.__wbg_ptr, level, col, row);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Tiles covering a pixel window on `level`, as a JSON array of
     * `{col,row,offset,length}`. Fetch each byte range, then `decode_tile_f64`.
     * @param {number} level
     * @param {number} x
     * @param {number} y
     * @param {number} w
     * @param {number} h
     * @returns {string}
     */
    tiles_for_window(level, x, y, w, h) {
        let deferred2_0;
        let deferred2_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.cogstream_tiles_for_window(retptr, this.__wbg_ptr, level, x, y, w, h);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            var ptr1 = r0;
            var len1 = r1;
            if (r3) {
                ptr1 = 0; len1 = 0;
                throw takeObject(r2);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred2_0, deferred2_1, 1);
        }
    }
}
if (Symbol.dispose) CogStream.prototype[Symbol.dispose] = CogStream.prototype.free;

/**
 * A parsed GeoTIFF held in memory. Construct once, then call the accessor and
 * `read_*` methods many times without re-parsing the file.
 */
export class GeoTiffReader {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        GeoTiffReaderFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_geotiffreader_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    get bands() {
        const ret = wasm.geotiffreader_bands(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    get bits_per_sample() {
        const ret = wasm.geotiffreader_bits_per_sample(this.__wbg_ptr);
        return ret;
    }
    /**
     * Bounding box as `[min_x, min_y, max_x, max_y]`, or empty if not georeferenced.
     * @returns {Float64Array}
     */
    bounding_box() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_bounding_box(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Bounds `[min_lon, min_lat, max_lon, max_lat]` in WGS84 degrees, or empty
     * if not convertible.
     * @returns {Float64Array}
     */
    bounds_lonlat() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_bounds_lonlat(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Image center `[x, y]` in the dataset CRS, or empty if not georeferenced.
     * @returns {Float64Array}
     */
    center() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_center(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Image center `[lon, lat]` in WGS84 degrees, or empty if not georeferenced
     * or the CRS is not convertible.
     * @returns {Float64Array}
     */
    center_lonlat() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_center_lonlat(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * @returns {string}
     */
    get compression() {
        let deferred1_0;
        let deferred1_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_compression(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            deferred1_0 = r0;
            deferred1_1 = r1;
            return getStringFromWasm0(r0, r1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * EPSG code, or `undefined` if the file is not georeferenced by EPSG.
     * @returns {number | undefined}
     */
    get epsg() {
        const ret = wasm.geotiffreader_epsg(this.__wbg_ptr);
        return ret === Number.MAX_SAFE_INTEGER ? undefined : ret;
    }
    /**
     * Affine geo-transform as `[x_origin, pixel_width, row_rotation,
     * y_origin, col_rotation, pixel_height]`, or an empty array if absent.
     * @returns {Float64Array}
     */
    geo_transform() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_geo_transform(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * @returns {number}
     */
    get height() {
        const ret = wasm.geotiffreader_height(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Full metadata as a JSON string (same shape as [`geotiff_info`]).
     * @returns {string}
     */
    info_json() {
        let deferred1_0;
        let deferred1_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_info_json(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            deferred1_0 = r0;
            deferred1_1 = r1;
            return getStringFromWasm0(r0, r1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @returns {boolean}
     */
    get is_bigtiff() {
        const ret = wasm.geotiffreader_is_bigtiff(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Parse a GeoTIFF / BigTIFF / COG from raw bytes.
     * @param {Uint8Array} data
     */
    constructor(data) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.geotiffreader_new(retptr, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            if (r2) {
                throw takeObject(r1);
            }
            this.__wbg_ptr = r0;
            GeoTiffReaderFinalization.register(this, this.__wbg_ptr, this);
            return this;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * No-data sentinel, or `undefined` if none is declared.
     * @returns {number | undefined}
     */
    get nodata() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_nodata(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r2 = getDataViewMemory0().getFloat64(retptr + 8 * 1, true);
            return r0 === 0 ? undefined : r2;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Read every band as `f64`, interleaved per pixel (`band0,band1,...`).
     * @returns {Float64Array}
     */
    read_all_f64() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_all_f64(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Read a band's raw, undecoded-to-native bytes. `Uint8Array`.
     * @param {number} band
     * @returns {Uint8Array}
     */
    read_band_bytes(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_bytes(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 1, 1);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Native `f32` band. `Float32Array`.
     * @param {number} band
     * @returns {Float32Array}
     */
    read_band_f32(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_f32(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayF32FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 4, 4);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Read a band as `f64`, converting from any on-disk type. `Float64Array`.
     * @param {number} band
     * @returns {Float64Array}
     */
    read_band_f64(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_f64(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Native `i16` band. `Int16Array`.
     * @param {number} band
     * @returns {Int16Array}
     */
    read_band_i16(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_i16(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayI16FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 2, 2);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Native `i32` band. `Int32Array`.
     * @param {number} band
     * @returns {Int32Array}
     */
    read_band_i32(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_i32(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayI32FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 4, 4);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Native `i8` band. `Int8Array`.
     * @param {number} band
     * @returns {Int8Array}
     */
    read_band_i8(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_i8(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayI8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 1, 1);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Native `u16` band. `Uint16Array`.
     * @param {number} band
     * @returns {Uint16Array}
     */
    read_band_u16(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_u16(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayU16FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 2, 2);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Native `u32` band. `Uint32Array`.
     * @param {number} band
     * @returns {Uint32Array}
     */
    read_band_u32(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_u32(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayU32FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 4, 4);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Native `u8` band. `Uint8Array`.
     * @param {number} band
     * @returns {Uint8Array}
     */
    read_band_u8(band) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_read_band_u8(retptr, this.__wbg_ptr, band);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 1, 1);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * @returns {string}
     */
    get sample_format() {
        let deferred1_0;
        let deferred1_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_sample_format(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            deferred1_0 = r0;
            deferred1_1 = r1;
            return getStringFromWasm0(r0, r1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Band-0 statistics as a JSON string (same shape as [`geotiff_stats`]).
     * @returns {string}
     */
    stats_json() {
        let deferred1_0;
        let deferred1_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_stats_json(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            deferred1_0 = r0;
            deferred1_1 = r1;
            return getStringFromWasm0(r0, r1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * GDAL value transform as `[scale, offset]` (physical = raw*scale+offset),
     * or empty if none. Apply to `read_*` outputs to get physical values.
     * @returns {Float64Array}
     */
    value_transform() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.geotiffreader_value_transform(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var v1 = getArrayF64FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 8, 8);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * @returns {number}
     */
    get width() {
        const ret = wasm.geotiffreader_width(this.__wbg_ptr);
        return ret >>> 0;
    }
}
if (Symbol.dispose) GeoTiffReader.prototype[Symbol.dispose] = GeoTiffReader.prototype.free;

export class PmtilesExtractor {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        PmtilesExtractorFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_pmtilesextractor_free(ptr, 0);
    }
    /**
     * True once every needed range has been fed; `finish()` is then valid.
     * @returns {boolean}
     */
    get done() {
        const ret = wasm.pmtilesextractor_done(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Hand back the bytes of one `wanted_json()` range, identified by its
     * offset. Ranges may be fed in any order.
     * @param {number} offset
     * @param {Uint8Array} bytes
     */
    feed(offset, bytes) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_export2);
            const len0 = WASM_VECTOR_LEN;
            wasm.pmtilesextractor_feed(retptr, this.__wbg_ptr, offset, ptr0, len0);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            if (r1) {
                throw takeObject(r0);
            }
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Assemble the extracted archive. Consumes the extractor's buffers; the
     * returned `Uint8Array` is a complete `.pmtiles` file.
     * @returns {Uint8Array}
     */
    finish() {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.pmtilesextractor_finish(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            if (r3) {
                throw takeObject(r2);
            }
            var v1 = getArrayU8FromWasm0(r0, r1).slice();
            wasm.__wbindgen_export(r0, r1 * 1, 1);
            return v1;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Source archive header as JSON (`{}` until the first feed): zooms,
     * bounds, tile type/compression, tile counts. Lets a UI validate the
     * request and describe the source before committing to the download.
     * @returns {string}
     */
    header_json() {
        let deferred2_0;
        let deferred2_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.pmtilesextractor_header_json(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            var ptr1 = r0;
            var len1 = r1;
            if (r3) {
                ptr1 = 0; len1 = 0;
                throw takeObject(r2);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * Plan an extraction of `min_zoom..=max_zoom` tiles intersecting the
     * WGS84 bbox. Zooms are clamped to what the source archive contains once
     * its header arrives; `min_zoom` 0 keeps the basemap usable zoomed out.
     * @param {number} min_lon
     * @param {number} min_lat
     * @param {number} max_lon
     * @param {number} max_lat
     * @param {number} min_zoom
     * @param {number} max_zoom
     */
    constructor(min_lon, min_lat, max_lon, max_lat, min_zoom, max_zoom) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.pmtilesextractor_new(retptr, min_lon, min_lat, max_lon, max_lat, min_zoom, max_zoom);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            if (r2) {
                throw takeObject(r1);
            }
            this.__wbg_ptr = r0;
            PmtilesExtractorFinalization.register(this, this.__wbg_ptr, this);
            return this;
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Progress as JSON: `{"phase":"header|directories|data|done",
     * "tiles_selected":n,"blobs_total":n,"data_bytes_total":n,
     * "data_bytes_received":n,"estimated_output_bytes":n}`.
     * @returns {string}
     */
    progress_json() {
        let deferred2_0;
        let deferred2_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.pmtilesextractor_progress_json(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            var ptr1 = r0;
            var len1 = r1;
            if (r3) {
                ptr1 = 0; len1 = 0;
                throw takeObject(r2);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * Coalesce tile-data requests whose byte gap is at most this (default
     * 65,536). Larger values trade overfetch for fewer HTTP round-trips.
     * @param {number} max_gap
     */
    set_max_range_gap(max_gap) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.pmtilesextractor_set_max_range_gap(retptr, this.__wbg_ptr, max_gap);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            if (r1) {
                throw takeObject(r0);
            }
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Cap on addressed tiles (default 2,000,000). Raise for huge desktop
     * extracts; lower to fail fast in memory-constrained embeds.
     * @param {number} max_tiles
     */
    set_max_tiles(max_tiles) {
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.pmtilesextractor_set_max_tiles(retptr, this.__wbg_ptr, max_tiles);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            if (r1) {
                throw takeObject(r0);
            }
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
        }
    }
    /**
     * Outstanding byte ranges the host should fetch, as a JSON array of
     * `{"offset":n,"length":n}`. Empty array when nothing is outstanding.
     * @returns {string}
     */
    wanted_json() {
        let deferred2_0;
        let deferred2_1;
        try {
            const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
            wasm.pmtilesextractor_wanted_json(retptr, this.__wbg_ptr);
            var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
            var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
            var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
            var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
            var ptr1 = r0;
            var len1 = r1;
            if (r3) {
                ptr1 = 0; len1 = 0;
                throw takeObject(r2);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_add_to_stack_pointer(16);
            wasm.__wbindgen_export(deferred2_0, deferred2_1, 1);
        }
    }
}
if (Symbol.dispose) PmtilesExtractor.prototype[Symbol.dispose] = PmtilesExtractor.prototype.free;

/**
 * Install a panic hook so Rust panics surface as readable `console.error`
 * messages instead of an opaque `RuntimeError: unreachable`.
 */
export function __start() {
    wasm.__start();
}

/**
 * Convex hull of a 2D point set. Input is `[x0,y0,x1,y1,...]`; output is the
 * hull ring as `[x0,y0,...]` (closed). Needs at least 3 points.
 * @param {Float64Array} points_xy
 * @returns {Float64Array}
 */
export function convex_hull(points_xy) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArrayF64ToWasm0(points_xy, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        wasm.convex_hull(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        if (r3) {
            throw takeObject(r2);
        }
        var v2 = getArrayF64FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export(r0, r1 * 8, 8);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Decode only a GeoTIFF's header and return its metadata as JSON. O(header)
 * memory, so it works on multi-gigabyte rasters that whole-image reads cannot
 * fit in WASM's 4 GiB address space.
 *
 * `{"ok":true,"width","height","bands","epsg"|null,"nodata"|null,
 *   "bits_per_sample","sample_format","compression","tiled","bigtiff"}`
 * @param {Uint8Array} data
 * @returns {string}
 */
export function geotiff_info(data) {
    let deferred2_0;
    let deferred2_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        wasm.geotiff_info(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        deferred2_0 = r0;
        deferred2_1 = r1;
        return getStringFromWasm0(r0, r1);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Read a single band of pixel values as an `f64` `Float64Array` (any on-disk
 * sample format is converted), row-major, length `width * height`.
 * @param {Uint8Array} data
 * @param {number} band
 * @returns {Float64Array}
 */
export function geotiff_read_band_f64(data, band) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        wasm.geotiff_read_band_f64(retptr, ptr0, len0, band);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        if (r3) {
            throw takeObject(r2);
        }
        var v2 = getArrayF64FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export(r0, r1 * 8, 8);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Decode a GeoTIFF and return band-0 summary statistics as JSON:
 * `{"ok":true,"width","height","bands","epsg","valid","min","max","mean"}`.
 * @param {Uint8Array} data
 * @returns {string}
 */
export function geotiff_stats(data) {
    let deferred2_0;
    let deferred2_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        wasm.geotiff_stats(retptr, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        deferred2_0 = r0;
        deferred2_1 = r1;
        return getStringFromWasm0(r0, r1);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred2_0, deferred2_1, 1);
    }
}

/**
 * LiDAR formats this build can read from memory.
 * @returns {string}
 */
export function lidar_formats() {
    let deferred1_0;
    let deferred1_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        wasm.lidar_formats(retptr);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        deferred1_0 = r0;
        deferred1_1 = r1;
        return getStringFromWasm0(r0, r1);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Read a LiDAR file's metadata as JSON. For LAS/LAZ this is header-only (count,
 * bounds, CRS, point format, COPC flag) and never decodes points:
 * `{"ok":true,"format","points","epsg"|null,"point_format"|null,
 *   "bounds":[min_x,min_y,min_z,max_x,max_y,max_z]|null,"copc":bool}`.
 * @param {Uint8Array} data
 * @param {string} format
 * @returns {string}
 */
export function lidar_info(data, format) {
    let deferred4_0;
    let deferred4_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(format, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
        const len1 = WASM_VECTOR_LEN;
        wasm.lidar_info(retptr, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        var ptr3 = r0;
        var len3 = r1;
        if (r3) {
            ptr3 = 0; len3 = 0;
            throw takeObject(r2);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred4_0, deferred4_1, 1);
    }
}

/**
 * Read per-point classification codes as a `Uint8Array` (length `point_count`).
 * @param {Uint8Array} data
 * @param {string} format
 * @returns {Uint8Array}
 */
export function lidar_read_classification(data, format) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(format, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
        const len1 = WASM_VECTOR_LEN;
        wasm.lidar_read_classification(retptr, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        if (r3) {
            throw takeObject(r2);
        }
        var v3 = getArrayU8FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export(r0, r1 * 1, 1);
        return v3;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Read per-point intensity as a `Uint16Array` (length `point_count`).
 * @param {Uint8Array} data
 * @param {string} format
 * @returns {Uint16Array}
 */
export function lidar_read_intensity(data, format) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(format, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
        const len1 = WASM_VECTOR_LEN;
        wasm.lidar_read_intensity(retptr, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        if (r3) {
            throw takeObject(r2);
        }
        var v3 = getArrayU16FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export(r0, r1 * 2, 2);
        return v3;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Read all point coordinates as an interleaved `Float64Array`
 * `[x0,y0,z0, x1,y1,z1, ...]` (length `3 * point_count`).
 *
 * Guarded against 32-bit memory blowup; very large clouds return a clean error
 * (read the header with `lidar_info`, or downsample on your side).
 * @param {Uint8Array} data
 * @param {string} format
 * @returns {Float64Array}
 */
export function lidar_read_xyz(data, format) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(format, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
        const len1 = WASM_VECTOR_LEN;
        wasm.lidar_read_xyz(retptr, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        if (r3) {
            throw takeObject(r2);
        }
        var v3 = getArrayF64FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export(r0, r1 * 8, 8);
        return v3;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Global Moran's I spatial autocorrelation for point data, using a binary
 * distance-band spatial weights matrix (neighbors within `distance_threshold`).
 *
 * `points_xy` is `[x0,y0,...]`, `values` is one value per point. Returns JSON:
 * `{"ok":true,"morans_i","expected","variance","z_score","p_value","n"}`.
 *
 * Builds neighbors in O(n^2); intended for up to a few thousand points.
 * @param {Float64Array} points_xy
 * @param {Float64Array} values
 * @param {number} distance_threshold
 * @returns {string}
 */
export function morans_i(points_xy, values, distance_threshold) {
    let deferred4_0;
    let deferred4_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArrayF64ToWasm0(points_xy, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(values, wasm.__wbindgen_export2);
        const len1 = WASM_VECTOR_LEN;
        wasm.morans_i(retptr, ptr0, len0, ptr1, len1, distance_threshold);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        var ptr3 = r0;
        var len3 = r1;
        if (r3) {
            ptr3 = 0; len3 = 0;
            throw takeObject(r2);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred4_0, deferred4_1, 1);
    }
}

/**
 * Reproject a bbox between two EPSG CRSs.
 *
 * Input and output order is `[min_x, min_y, max_x, max_y]`. The bbox edges are
 * densified so projected extrema that fall along an edge are preserved better
 * than a corner-only transform.
 * @param {number} src_epsg
 * @param {number} dst_epsg
 * @param {Float64Array} bbox
 * @returns {Float64Array}
 */
export function transform_bbox_epsg(src_epsg, dst_epsg, bbox) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArrayF64ToWasm0(bbox, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        wasm.transform_bbox_epsg(retptr, src_epsg, dst_epsg, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        if (r3) {
            throw takeObject(r2);
        }
        var v2 = getArrayF64FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export(r0, r1 * 8, 8);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Reproject x,y coordinate pairs between two EPSG CRSs.
 * @param {number} src_epsg
 * @param {number} dst_epsg
 * @param {Float64Array} xy
 * @returns {Float64Array}
 */
export function transform_points_epsg(src_epsg, dst_epsg, xy) {
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArrayF64ToWasm0(xy, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        wasm.transform_points_epsg(retptr, src_epsg, dst_epsg, ptr0, len0);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        if (r3) {
            throw takeObject(r2);
        }
        var v2 = getArrayF64FromWasm0(r0, r1).slice();
        wasm.__wbindgen_export(r0, r1 * 8, 8);
        return v2;
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
    }
}

/**
 * Vector formats this build can read from memory (comma-separated).
 * @returns {string}
 */
export function vector_formats() {
    let deferred1_0;
    let deferred1_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        wasm.vector_formats(retptr);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        deferred1_0 = r0;
        deferred1_1 = r1;
        return getStringFromWasm0(r0, r1);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Read a vector dataset and return metadata as JSON:
 * `{"ok":true,"name","features","geometry","epsg"|null,"fields":[...],
 *   "bbox":[min_x,min_y,max_x,max_y]|null}`.
 * @param {Uint8Array} data
 * @param {string} format
 * @returns {string}
 */
export function vector_info(data, format) {
    let deferred4_0;
    let deferred4_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(format, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
        const len1 = WASM_VECTOR_LEN;
        wasm.vector_info(retptr, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        var ptr3 = r0;
        var len3 = r1;
        if (r3) {
            ptr3 = 0; len3 = 0;
            throw takeObject(r2);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred4_0, deferred4_1, 1);
    }
}

/**
 * Read a vector dataset and return it as a GeoJSON `FeatureCollection` string.
 * @param {Uint8Array} data
 * @param {string} format
 * @returns {string}
 */
export function vector_to_geojson(data, format) {
    let deferred4_0;
    let deferred4_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(format, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
        const len1 = WASM_VECTOR_LEN;
        wasm.vector_to_geojson(retptr, ptr0, len0, ptr1, len1);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        var ptr3 = r0;
        var len3 = r1;
        if (r3) {
            ptr3 = 0; len3 = 0;
            throw takeObject(r2);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred4_0, deferred4_1, 1);
    }
}

/**
 * Read a vector dataset, reproject it to `dst_epsg`, and return GeoJSON.
 * Uses the bundled pure-Rust projection engine (full EPSG support).
 *
 * `src_epsg` overrides the source CRS: pass `0` to use the layer's own CRS, or
 * fall back to EPSG:4326 if it declares none (GeoJSON is WGS84 by RFC 7946).
 * @param {Uint8Array} data
 * @param {string} format
 * @param {number} dst_epsg
 * @param {number} src_epsg
 * @returns {string}
 */
export function vector_to_geojson_reproject(data, format, dst_epsg, src_epsg) {
    let deferred4_0;
    let deferred4_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        const ptr0 = passArray8ToWasm0(data, wasm.__wbindgen_export2);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(format, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
        const len1 = WASM_VECTOR_LEN;
        wasm.vector_to_geojson_reproject(retptr, ptr0, len0, ptr1, len1, dst_epsg, src_epsg);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        var r2 = getDataViewMemory0().getInt32(retptr + 4 * 2, true);
        var r3 = getDataViewMemory0().getInt32(retptr + 4 * 3, true);
        var ptr3 = r0;
        var len3 = r1;
        if (r3) {
            ptr3 = 0; len3 = 0;
            throw takeObject(r2);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred4_0, deferred4_1, 1);
    }
}

/**
 * Semantic version of this crate, exposed for runtime feature detection.
 * @returns {string}
 */
export function version() {
    let deferred1_0;
    let deferred1_1;
    try {
        const retptr = wasm.__wbindgen_add_to_stack_pointer(-16);
        wasm.version(retptr);
        var r0 = getDataViewMemory0().getInt32(retptr + 4 * 0, true);
        var r1 = getDataViewMemory0().getInt32(retptr + 4 * 1, true);
        deferred1_0 = r0;
        deferred1_1 = r1;
        return getStringFromWasm0(r0, r1);
    } finally {
        wasm.__wbindgen_add_to_stack_pointer(16);
        wasm.__wbindgen_export(deferred1_0, deferred1_1, 1);
    }
}
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_throw_ea4887a5f8f9a9db: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbg_error_a6fa202b58aa1cd3: function(arg0, arg1) {
            let deferred0_0;
            let deferred0_1;
            try {
                deferred0_0 = arg0;
                deferred0_1 = arg1;
                console.error(getStringFromWasm0(arg0, arg1));
            } finally {
                wasm.__wbindgen_export(deferred0_0, deferred0_1, 1);
            }
        },
        __wbg_new_227d7c05414eb861: function() {
            const ret = new Error();
            return addHeapObject(ret);
        },
        __wbg_stack_3b0d974bbf31e44f: function(arg0, arg1) {
            const ret = getObject(arg1).stack;
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_export2, wasm.__wbindgen_export3);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbindgen_cast_0000000000000001: function(arg0, arg1) {
            // Cast intrinsic for `Ref(String) -> Externref`.
            const ret = getStringFromWasm0(arg0, arg1);
            return addHeapObject(ret);
        },
        __wbindgen_object_drop_ref: function(arg0) {
            takeObject(arg0);
        },
    };
    return {
        __proto__: null,
        "./geolibre_wasm_bg.js": import0,
    };
}

const CogBuilderFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_cogbuilder_free(ptr, 1));
const CogStreamFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_cogstream_free(ptr, 1));
const GeoTiffReaderFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_geotiffreader_free(ptr, 1));
const PmtilesExtractorFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_pmtilesextractor_free(ptr, 1));

function addHeapObject(obj) {
    if (heap_next === heap.length) heap.push(heap.length + 1);
    const idx = heap_next;
    heap_next = heap[idx];

    heap[idx] = obj;
    return idx;
}

function dropObject(idx) {
    if (idx < 1028) return;
    heap[idx] = heap_next;
    heap_next = idx;
}

function getArrayF32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}

function getArrayI16FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getInt16ArrayMemory0().subarray(ptr / 2, ptr / 2 + len);
}

function getArrayI32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getInt32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayI8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getInt8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

function getArrayU16FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint16ArrayMemory0().subarray(ptr / 2, ptr / 2 + len);
}

function getArrayU32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

let cachedDataViewMemory0 = null;
function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

let cachedFloat32ArrayMemory0 = null;
function getFloat32ArrayMemory0() {
    if (cachedFloat32ArrayMemory0 === null || cachedFloat32ArrayMemory0.byteLength === 0) {
        cachedFloat32ArrayMemory0 = new Float32Array(wasm.memory.buffer);
    }
    return cachedFloat32ArrayMemory0;
}

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
}

let cachedInt16ArrayMemory0 = null;
function getInt16ArrayMemory0() {
    if (cachedInt16ArrayMemory0 === null || cachedInt16ArrayMemory0.byteLength === 0) {
        cachedInt16ArrayMemory0 = new Int16Array(wasm.memory.buffer);
    }
    return cachedInt16ArrayMemory0;
}

let cachedInt32ArrayMemory0 = null;
function getInt32ArrayMemory0() {
    if (cachedInt32ArrayMemory0 === null || cachedInt32ArrayMemory0.byteLength === 0) {
        cachedInt32ArrayMemory0 = new Int32Array(wasm.memory.buffer);
    }
    return cachedInt32ArrayMemory0;
}

let cachedInt8ArrayMemory0 = null;
function getInt8ArrayMemory0() {
    if (cachedInt8ArrayMemory0 === null || cachedInt8ArrayMemory0.byteLength === 0) {
        cachedInt8ArrayMemory0 = new Int8Array(wasm.memory.buffer);
    }
    return cachedInt8ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    return decodeText(ptr >>> 0, len);
}

let cachedUint16ArrayMemory0 = null;
function getUint16ArrayMemory0() {
    if (cachedUint16ArrayMemory0 === null || cachedUint16ArrayMemory0.byteLength === 0) {
        cachedUint16ArrayMemory0 = new Uint16Array(wasm.memory.buffer);
    }
    return cachedUint16ArrayMemory0;
}

let cachedUint32ArrayMemory0 = null;
function getUint32ArrayMemory0() {
    if (cachedUint32ArrayMemory0 === null || cachedUint32ArrayMemory0.byteLength === 0) {
        cachedUint32ArrayMemory0 = new Uint32Array(wasm.memory.buffer);
    }
    return cachedUint32ArrayMemory0;
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function getObject(idx) { return heap[idx]; }

let heap = new Array(1024).fill(undefined);
heap.push(undefined, null, true, false);

let heap_next = heap.length;

function passArray32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getUint32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passArray8ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 1, 1) >>> 0;
    getUint8ArrayMemory0().set(arg, ptr / 1);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passArrayF32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getFloat32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passArrayF64ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 8, 8) >>> 0;
    getFloat64ArrayMemory0().set(arg, ptr / 8);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeObject(idx) {
    const ret = getObject(idx);
    dropObject(idx);
    return ret;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();

if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function (arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {
            read: arg.length,
            written: buf.length
        };
    };
}

let WASM_VECTOR_LEN = 0;

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    cachedDataViewMemory0 = null;
    cachedFloat32ArrayMemory0 = null;
    cachedFloat64ArrayMemory0 = null;
    cachedInt16ArrayMemory0 = null;
    cachedInt32ArrayMemory0 = null;
    cachedInt8ArrayMemory0 = null;
    cachedUint16ArrayMemory0 = null;
    cachedUint32ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = module.ok && expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('geolibre_wasm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
