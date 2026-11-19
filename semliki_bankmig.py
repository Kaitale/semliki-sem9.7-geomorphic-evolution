"""semliki_bankmig.py -- bank retreat and centreline migration rates from annual channel-water masks.

Part of the reproducibility package of Nafula, M. and Kakooza, C. (2026), "Meander migration, cut-off risk and
bed stability on the lower Semliki River, Albertine Rift", preprint, https://doi.org/10.5281/zenodo.24367442.

Method (DSAS-style transects on a fixed baseline, applied to a river):
  1  each year's water mask is reduced to the CHANNEL: the connected water body that touches the baseline
     (oxbow lakes, ponds and swamp pools drop out; islands smaller than fill_ha are filled);
  2  the baseline is resampled every `spacing` m and smoothed; a transect normal to it is cast at every
     station, +-`half` m;
  3  on every transect, for every year, the left and right bank positions are found at sub-pixel precision
     (the 0.5 crossing of the bilinearly sampled mask) on the water run that continues last year's channel;
     years whose transect crosses cloud/no-data ("valid" band = 0) are dropped for that transect;
  4  per transect and bank: linear-regression rate (m/yr) with its 95 % confidence interval, end-point rate,
     net movement; retreat is positive when the bank moves away from the channel (erosion);
  5  bends are the reaches between curvature inflections; per bend the OUTER bank's rates are summarised,
     with an area-based check (eroded area / bend length / years);
  6  meander necks closer than neck_factor x channel width across land, but far apart along the river,
     are listed as cut-off candidates;
  7  setback per transect = upper-95 % retreat rate x design life + 2 x positional uncertainty +
     bank-stability allowance (bank height / tan(stable angle)).

Positional uncertainty of a bank in one year: e = sqrt((pixel/2)^2 + registration^2).
Licence: MIT (see LICENSE).
"""
import argparse, glob, json, math, os, re, sys
import numpy as np
from scipy import ndimage, stats
from scipy.spatial import cKDTree



# ------------------------------------------------------------------ raster IO (rasterio or GDAL)
def read_mask(path):
    """returns water (uint8), valid (uint8), transform (a, b, c, d, e, f) with x = a*col + b*row + c, y = d*col + e*row + f.
    .npz files (gee_water_masks.py --download-dir) need nothing beyond numpy; GeoTIFFs need rasterio or GDAL."""
    if path.lower().endswith(".npz"):
        z = np.load(path)
        return z["water"].astype(np.uint8), z["valid"].astype(np.uint8), tuple(float(v) for v in z["transform"])
    try:
        import rasterio
        with rasterio.open(path) as r:
            w = r.read(1).astype(np.uint8)
            v = r.read(2).astype(np.uint8) if r.count > 1 else np.ones_like(w)
            t = r.transform
            return w, v, (t.a, t.b, t.c, t.d, t.e, t.f)
    except ImportError:
        from osgeo import gdal
        ds = gdal.Open(path)
        w = ds.GetRasterBand(1).ReadAsArray().astype(np.uint8)
        v = ds.GetRasterBand(2).ReadAsArray().astype(np.uint8) if ds.RasterCount > 1 else np.ones_like(w)
        c, a, b, f, d, e = ds.GetGeoTransform()
        return w, v, (a, b, c, d, e, f)


def year_of(path):
    m = re.findall(r"(19\d\d|20\d\d)", os.path.basename(path))
    return int(m[-1])


# ------------------------------------------------------------------ geometry helpers
def to_pix(T, x, y):
    a, b, c, d, e, f = T
    return (y - f) / e, (x - c) / a          # row, col (north-up rasters)


def resample_line(xy, spacing):
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    s = np.arange(0, seg[-1], spacing)
    return np.column_stack([np.interp(s, seg, xy[:, 0]), np.interp(s, seg, xy[:, 1])]), s


def smooth_line(xy, window_pts):
    if window_pts < 3:
        return xy
    k = np.ones(window_pts) / window_pts
    pad = window_pts // 2
    out = np.empty_like(xy)
    for j in range(2):
        p = np.r_[np.full(pad, xy[0, j]), xy[:, j], np.full(pad, xy[-1, j])]
        out[:, j] = np.convolve(p, k, mode="valid")[:len(xy)]
    out[0], out[-1] = xy[0], xy[-1]
    return out


def normals(xy):
    t = np.gradient(xy, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    return np.column_stack([t[:, 1], -t[:, 0]])          # right-hand normal (to the right looking downstream)


def curvature(xy, s):
    dx, dy = np.gradient(xy[:, 0], s), np.gradient(xy[:, 1], s)
    ddx, ddy = np.gradient(dx, s), np.gradient(dy, s)
    return (dx * ddy - dy * ddx) / np.power(dx * dx + dy * dy, 1.5)   # + = turning left


def baseline_from_masks(channel_stack, T):
    """skeleton of the multi-year channel union, longest geodesic path; returned upstream?->downstream unknown:
    the caller must check the direction (pass --reverse if the first point is at the downstream end)"""
    from skimage.morphology import skeletonize
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra
    u = ndimage.binary_closing(channel_stack.any(axis=0), iterations=2)
    sk = skeletonize(u)
    rr, cc = np.nonzero(sk)
    idx = -np.ones(sk.shape, int)
    idx[rr, cc] = np.arange(len(rr))
    rows, cols, wts = [], [], []
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == dc == 0:
                continue
            r2, c2 = rr + dr, cc + dc
            ok = (r2 >= 0) & (r2 < sk.shape[0]) & (c2 >= 0) & (c2 < sk.shape[1])
            j = np.full(len(rr), -1)
            j[ok] = idx[r2[ok], c2[ok]]
            m = j >= 0
            rows += list(np.nonzero(m)[0]); cols += list(j[m]); wts += [math.hypot(dr, dc)] * int(m.sum())
    G = csr_matrix((wts, (rows, cols)), shape=(len(rr), len(rr)))
    d0 = dijkstra(G, indices=0)
    a = int(np.nanargmax(np.where(np.isfinite(d0), d0, -1)))
    da, pred = dijkstra(G, indices=a, return_predecessors=True)
    b = int(np.nanargmax(np.where(np.isfinite(da), da, -1)))
    path = [b]
    while path[-1] != a:
        path.append(pred[path[-1]])
    a_, b_, c_, d_, e_, f_ = T
    pr, pc = rr[path] + 0.5, cc[path] + 0.5
    return np.column_stack([a_ * pc + b_ * pr + c_, d_ * pc + e_ * pr + f_])


# ------------------------------------------------------------------ channel cleaning
def channel_only(water, T, baseline_xy, fill_ha, pixel):
    w = ndimage.binary_closing(water.astype(bool), structure=np.ones((3, 3)), iterations=1)
    holes = ndimage.binary_fill_holes(w) & ~w
    lab, n = ndimage.label(holes)
    if n:
        sizes = ndimage.sum(holes, lab, range(1, n + 1)) * pixel * pixel / 1e4
        small = np.isin(lab, np.nonzero(sizes < fill_ha)[0] + 1)
        w |= small
    lab, n = ndimage.label(w, structure=np.ones((3, 3)))
    r, c = to_pix(T, baseline_xy[:, 0], baseline_xy[:, 1])
    r, c = np.round(r).astype(int), np.round(c).astype(int)
    ok = (r >= 0) & (r < w.shape[0]) & (c >= 0) & (c < w.shape[1])
    # components on the baseline or within 3 pixels of it (a channel that has drifted off the baseline)
    near = np.zeros_like(w)
    near[r[ok], c[ok]] = True
    near = ndimage.binary_dilation(near, iterations=3)
    seeds = set(np.unique(lab[near & w])) - {0}
    return np.isin(lab, list(seeds))


# ------------------------------------------------------------------ banks on transects
def bank_positions(chan_f, valid, T, pts, nrm, half, step, prev, gap_px, pixel):
    """for each transect: (left, right) bank offsets in m (negative = left of baseline), NaN if not found/invalid.
    prev: previous year's centre offsets (to follow the same channel)"""
    s = np.arange(-half, half + step, step)
    X = pts[:, 0][:, None] + nrm[:, 0][:, None] * s[None, :]
    Y = pts[:, 1][:, None] + nrm[:, 1][:, None] * s[None, :]
    r, c = to_pix(T, X, Y)
    prof = ndimage.map_coordinates(chan_f, [r - 0.5, c - 0.5], order=1, mode="constant", cval=0.0)
    vprof = ndimage.map_coordinates(valid, [r - 0.5, c - 0.5], order=0, mode="constant", cval=0)
    out = np.full((len(pts), 2), np.nan)
    gap = max(1, int(round(gap_px * pixel / step)))
    for i in range(len(pts)):
        wet = prof[i] >= 0.5
        if not wet.any():
            continue
        # merge runs separated by short dry gaps (mid-channel bars)
        wet = ndimage.binary_closing(wet, structure=np.ones(2 * gap + 1)) | wet
        lab, n = ndimage.label(wet)
        target = prev[i] if np.isfinite(prev[i]) else 0.0
        best, bd = None, np.inf
        for k in range(1, n + 1):
            j = np.nonzero(lab == k)[0]
            d = 0.0 if s[j[0]] <= target <= s[j[-1]] else min(abs(s[j[0]] - target), abs(s[j[-1]] - target))
            if d < bd:
                best, bd = j, d
        j0, j1 = best[0], best[-1]
        if j0 == 0 or j1 == len(s) - 1:
            continue                                    # run reaches the transect end: half-length too short
        # sub-pixel crossings of 0.5
        def cross(a, b):
            pa, pb = prof[i][a], prof[i][b]
            return s[a] + (0.5 - pa) / (pb - pa) * (s[b] - s[a]) if pb != pa else s[a]
        left, right = cross(j0 - 1, j0), cross(j1, j1 + 1)
        if vprof[i][max(j0 - 3, 0):j1 + 4].min() == 0:
            continue                                    # cloud or no data on this transect this year
        out[i] = (left, right)
    return out


# ------------------------------------------------------------------ rates
def lrr(years, pos):
    m = np.isfinite(pos)
    if m.sum() < 3:
        return dict(n=int(m.sum()), rate=np.nan, ci95=np.nan, r2=np.nan, epr=np.nan, nsm=np.nan)
    y, p = years[m].astype(float), pos[m]
    res = stats.linregress(y, p)
    tcrit = stats.t.ppf(0.975, m.sum() - 2)
    return dict(n=int(m.sum()), rate=res.slope, ci95=tcrit * res.stderr, r2=res.rvalue ** 2,
                epr=(p[-1] - p[0]) / (y[-1] - y[0]), nsm=p[-1] - p[0])


def analyse(masks, T, baseline_xy, spacing=25.0, half=600.0, fill_ha=0.5, gap_px=2, smooth_m=200.0,
            registration=10.0, design_life=50.0, bank_height=4.0, stable_angle=26.6, neck_factor=3.0, pixel=None,
            presmoothed=False):
    """masks: {year: (water, valid)} on one grid.  Returns dict of tables (lists of dicts) and arrays."""
    years = np.array(sorted(masks))
    pixel = pixel or abs(T[0])
    step = pixel / 4.0
    shape = next(iter(masks.values()))[0].shape
    bxy = np.asarray(baseline_xy, float)
    r_, c_ = to_pix(T, bxy[:, 0], bxy[:, 1])
    inside = (r_ >= half / pixel) & (r_ < shape[0] - half / pixel) & (c_ >= half / pixel) & (c_ < shape[1] - half / pixel)
    if inside.sum() < 2:
        raise ValueError("the baseline does not cross the mask extent")
    k0, k1 = np.nonzero(inside)[0][[0, -1]]
    if presmoothed:       # the baseline is already resampled and smoothed (e.g. the published transects.csv x, y)
        base = bxy[k0:k1 + 1]
        s_base = np.r_[0, np.cumsum(np.hypot(*np.diff(base, axis=0).T))]
    else:
        base, s_base = resample_line(bxy[k0:k1 + 1], spacing)      # baseline clipped to the rasters
        base = smooth_line(base, max(3, int(smooth_m / spacing) | 1))
    nrm = normals(base)
    kap = curvature(smooth_line(base, 5), s_base)
    pos = np.full((len(years), len(base), 2), np.nan)
    chans = {}
    prev = np.zeros(len(base))
    for k, y in enumerate(years):
        water, valid = masks[y]
        ch = channel_only(water, T, base, fill_ha, pixel)
        chans[y] = ch
        bp = bank_positions(ch.astype(np.float32), valid, T, base, nrm, half, step, prev, gap_px, pixel)
        pos[k] = bp
        c = bp.mean(axis=1)
        prev = np.where(np.isfinite(c), c, prev)
    e_pos = math.sqrt((pixel / 2) ** 2 + registration ** 2)
    allow = bank_height / math.tan(math.radians(stable_angle))
    rows = []
    for i in range(len(base)):
        L, R = pos[:, i, 0], pos[:, i, 1]
        rl, rr_, rc = lrr(years, -L), lrr(years, R), lrr(years, (L + R) / 2)     # retreat positive = erosion
        w = np.nanmedian(R - L) if np.isfinite(R - L).any() else np.nan
        row = dict(station_m=round(float(s_base[i]), 1), x=round(float(base[i, 0]), 1), y=round(float(base[i, 1]), 1),
                   curvature_per_km=round(float(kap[i]) * 1000, 3), outer_bank=("right" if kap[i] > 0 else "left"),
                   width_m=round(float(w), 1) if np.isfinite(w) else None)
        for lab, r in (("left", rl), ("right", rr_), ("centre", rc)):
            row.update({f"{lab}_n": r["n"], f"{lab}_rate": _r(r["rate"]), f"{lab}_ci95": _r(r["ci95"]), f"{lab}_r2": _r(r["r2"], 2),
                        f"{lab}_epr": _r(r["epr"]), f"{lab}_nsm": _r(r["nsm"])})
            if lab != "centre":
                sig = np.isfinite(r["rate"]) and (r["rate"] - r["ci95"] > 0)
                row[f"{lab}_eroding"] = bool(sig)
                up = (r["rate"] + r["ci95"]) if np.isfinite(r["rate"]) else np.nan
                row[f"{lab}_setback_m"] = _r(max(up, 0.0) * design_life + 2 * e_pos + allow, 0) if np.isfinite(up) else None
        rows.append(row)
    bends = bends_table(rows, kap, s_base, years, chans, base, nrm, T, pixel)
    necks = necks_table(base, s_base, rows, neck_factor)
    return dict(years=years, base=base, normals=nrm, s=s_base, pos=pos, transects=rows, bends=bends, necks=necks,
                e_pos=e_pos, allowance=allow, chans=chans)


def _r(v, nd=2):
    return None if v is None or not np.isfinite(v) else round(float(v), nd)


def bends_table(rows, kap, s, years, chans, base, nrm, T, pixel, min_len=300.0):
    sign = np.sign(kap)
    cuts = [0] + [i for i in range(1, len(sign)) if sign[i] != sign[i - 1]] + [len(sign)]
    out = []
    y0, y1 = years[0], years[-1]
    ero = chans[y1] & ~chans[y0]
    for b, (i0, i1) in enumerate(zip(cuts[:-1], cuts[1:])):
        if s[min(i1, len(s) - 1)] - s[i0] < min_len:
            continue
        seg = rows[i0:i1]
        side = "right" if np.mean(kap[i0:i1]) > 0 else "left"
        rates = np.array([r[f"{side}_rate"] if r[f"{side}_rate"] is not None else np.nan for r in seg])
        sb = np.array([r[f"{side}_setback_m"] if r[f"{side}_setback_m"] is not None else np.nan for r in seg])
        apex = i0 + int(np.argmax(np.abs(kap[i0:i1])))
        # area-based check: eroded pixels inside this bend's transect fan, on the outer side
        a, b_, c, d, e, f = T
        rr, cc = np.nonzero(ero)
        X, Y = a * (cc + 0.5) + c, e * (rr + 0.5) + f
        tree = cKDTree(base)
        _, j = tree.query(np.column_stack([X, Y]))
        sd = np.sign(((X - base[j, 0]) * nrm[j, 0] + (Y - base[j, 1]) * nrm[j, 1]))
        inb = (j >= i0) & (j < i1) & (sd == (1 if side == "right" else -1))
        area = inb.sum() * pixel * pixel
        length = s[min(i1, len(s) - 1)] - s[i0]
        out.append(dict(bend=b + 1, partial=(i0 == 0 or i1 == len(s)), start_m=round(float(s[i0])), end_m=round(float(s[min(i1, len(s) - 1)])), length_m=round(float(length)),
                        apex_x=round(float(base[apex, 0])), apex_y=round(float(base[apex, 1])),
                        radius_m=round(1 / abs(kap[apex])) if kap[apex] != 0 else None, outer_bank=side,
                        n_transects=int(np.isfinite(rates).sum()),
                        outer_rate_median=_r(np.nanmedian(rates)) if np.isfinite(rates).any() else None,
                        outer_rate_p90=_r(np.nanpercentile(rates, 90)) if np.isfinite(rates).any() else None,
                        outer_rate_max=_r(np.nanmax(rates)) if np.isfinite(rates).any() else None,
                        area_based_rate=_r(area / length / (y1 - y0)),          # end-point estimate: more exposed to registration error than the LRR
                        setback_p90_m=_r(np.nanpercentile(sb, 90), 0) if np.isfinite(sb).any() else None))
    return out


def necks_table(base, s, rows, factor):
    w = np.nanmedian([r["width_m"] for r in rows if r["width_m"]]) if any(r["width_m"] for r in rows) else 100.0
    tree = cKDTree(base)
    pairs = tree.query_pairs(factor * w)
    out = []
    for i, j in pairs:
        i, j = min(i, j), max(i, j)
        along, across = s[j] - s[i], float(np.hypot(*(base[j] - base[i])))
        if along > 6 * factor * w and along > 5 * across:
            out.append((across, along, i, j))
    out.sort()
    keep, used = [], []
    for across, along, i, j in out:                      # one entry per neck (non-overlapping)
        if any(abs(i - a) < 20 and abs(j - b) < 20 for a, b in used):
            continue
        used.append((i, j))
        keep.append(dict(upstream_station_m=round(float(s[i])), downstream_station_m=round(float(s[j])),
                         neck_width_m=round(across - w), bypassed_length_m=round(float(along)),
                         x=round(float((base[i, 0] + base[j, 0]) / 2)), y=round(float((base[i, 1] + base[j, 1]) / 2))))
    return keep


# ------------------------------------------------------------------ outputs
def write_outputs(R, out, epsg=32635):
    import csv
    os.makedirs(out, exist_ok=True)
    for name in ("transects", "bends", "necks"):
        rows = R[name]
        with open(os.path.join(out, name + ".csv"), "w", newline="") as fh:
            if rows:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader(); w.writerows(rows)
            else:
                fh.write("none\n")
    with open(os.path.join(out, "bank_positions.csv"), "w") as fh:
        fh.write("station_m,year,left_m,right_m\n")
        for k, y in enumerate(R["years"]):
            for i, s in enumerate(R["s"]):
                L, Rr = R["pos"][k, i]
                if np.isfinite(L):
                    fh.write("%.1f,%d,%.2f,%.2f\n" % (s, y, L, Rr))
    feats = []
    for r, p, n in zip(R["transects"], R["base"], R["normals"]):
        feats.append({"type": "Feature", "properties": r,
                      "geometry": {"type": "LineString", "coordinates": [list(p - 300 * n), list(p + 300 * n)]}})
    json.dump({"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "EPSG:%d" % epsg}}, "features": feats},
              open(os.path.join(out, "transects.geojson"), "w"))
    try:
        import matplotlib.pyplot as plt
        s = R["s"] / 1000
        fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
        for lab, col in (("left", "#1f77b4"), ("right", "#d62828")):
            v = np.array([r[f"{lab}_rate"] if r[f"{lab}_rate"] is not None else np.nan for r in R["transects"]])
            ci = np.array([r[f"{lab}_ci95"] if r[f"{lab}_ci95"] is not None else np.nan for r in R["transects"]])
            ax[0].plot(s, v, color=col, lw=0.8, label=f"{lab} bank retreat (m/yr, + = erosion)")
            ax[0].fill_between(s, v - ci, v + ci, color=col, alpha=0.15)
        ax[0].axhline(0, color="#999", lw=0.6); ax[0].set_ylabel("m/yr"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.2)
        ax[1].plot(s, [r["curvature_per_km"] for r in R["transects"]], color="#444", lw=0.8)
        ax[1].set_ylabel("curvature (1/km, + = left turn)"); ax[1].set_xlabel("distance along baseline (km)"); ax[1].grid(alpha=0.2)
        ax[0].set_title("Bank retreat rates %d-%d (linear regression, 95%% band)" % (R["years"][0], R["years"][-1]), loc="left")
        fig.tight_layout(); fig.savefig(os.path.join(out, "fig_rates_along_river.png"), dpi=130); plt.close(fig)
        b = [x for x in R["bends"] if x["outer_rate_median"] is not None]
        if b:
            fig, ax = plt.subplots(figsize=(11, 3.6))
            xb = np.arange(len(b))
            ax.bar(xb, [x["outer_rate_p90"] for x in b], color="#d62828", alpha=0.35, label="outer bank, 90th percentile")
            ax.bar(xb, [x["outer_rate_median"] for x in b], color="#d62828", label="outer bank, median")
            ax.plot(xb, [x["area_based_rate"] for x in b], "k.", label="area-based check")
            ax.set_xticks(xb); ax.set_xticklabels([str(x["bend"]) for x in b], fontsize=7); ax.set_xlabel("bend (upstream to downstream)")
            ax.set_ylabel("retreat, m/yr"); ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.2)
            fig.tight_layout(); fig.savefig(os.path.join(out, "fig_bends.png"), dpi=130); plt.close(fig)
    except Exception as ex:                               # figures are a convenience
        print("figures skipped:", ex)


