#!/usr/bin/env python
"""
Extract and plot AirMSPI L1C multi-angle polarimetric data for one grid pixel.

AirMSPI L1C (GRP_TERRAIN) delivers one HDF5 file per viewing angle
(e.g. 646F, 467F, 000N, 467A, 646A). All views of one target are mapped onto
the same UTM grid, so a (x, y) grid index refers to the same ground location
in every file.

File layout used here
---------------------
Channel_Information/Center_wavelength            per channel (14)
Channel_Information/Solar_irradiance_at_1_AU     F0 per channel, W m-2 nm-1
HDFEOS/GRIDS/<band>nm_band/Data Fields/I         radiance, W m-2 sr-1 nm-1
HDFEOS/GRIDS/<band>nm_band/Data Fields/Q_scatter Q in scattering-plane frame
HDFEOS/GRIDS/<band>nm_band/Data Fields/U_scatter U in scattering-plane frame
HDFEOS/GRIDS/<band>nm_band/Data Fields/DOLP, Scattering_angle, View_zenith,
    View_azimuth, Sun_zenith, Sun_azimuth, Glint_angle, RDQI, I.mask, ...
HDFEOS/GRIDS/Ancillary/Data Fields/Latitude, Longitude, Elevation
HDFEOS/GRIDS/XDim, YDim                          UTM easting / northing (m)

Arrays are stored as [YDim, XDim], i.e. [row, column]. In this script
``x`` is the column index (easting) and ``y`` is the row index (northing,
increasing southward). The masks (I.mask, Q.mask, U.mask) are 1 where the
data are valid and 0 where there are no data.

Polarization (Q, U, DoLP) is available in the 470, 660 and 865 nm bands only;
I is available in all eight bands.

Normalization follows the "BRF equation" stored in the file attributes:
    BRF = radiance * pi * sun_distance^2 / (F0 * cos(sun_zenith))
and is applied the same way to Q and U. DoLP = sqrt(Q^2 + U^2) / I is
computed from the (box-averaged) Stokes parameters; for --box 1 it equals
the file's DOLP field.

Signed view zenith angle: positive for forward-looking views, negative for
aft-looking views. The sign is taken from the view azimuth relative to the
aircraft heading (azimuths in these files follow the photon-propagation
convention, which is consistent with the file's Scattering_angle field), with
the gimbal angle as fallback.

Examples
--------
    python airmspi_pixel_plot.py --x 2048 --y 2048
    python airmspi_pixel_plot.py --x 2048 --y 2048 --box 5 --bands pol
    python airmspi_pixel_plot.py --utm 330000 3845000 --csv pixel.csv
    python airmspi_pixel_plot.py --latlon 34.733 -112.853 --show
    python airmspi_pixel_plot.py --x 2048 --y 2048 --quicklook
"""

import argparse
import glob
import os
import re
import sys
from datetime import datetime, timedelta

import h5py
import numpy as np
import pandas as pd

DEFAULT_DATA_DIR = "/stor/z101/Data/FIREX-AQ/AirMSPI/L1C/2019/08/17"
DEFAULT_PATTERN = "AirMSPI_*.hdf"
GRID_ROOT = "HDFEOS/GRIDS"
FILE_ATTRS = "HDFEOS/ADDITIONAL/FILE_ATTRIBUTES"
FILL_VALUE = -999.0

# Plot styling. Each band keeps the same color and marker regardless of which
# bands are selected. The legend lists the polarimetric bands first, then the
# intensity-only bands; the colors were checked for color-vision-deficiency
# separation in that order.
BAND_STYLE = {
    "470": ("#2a78d6", "o"),  # blue
    "660": ("#eb6834", "s"),  # orange
    "865": ("#1baf7a", "D"),  # aqua
    "355": ("#6250d6", "^"),  # violet
    "380": ("#e87ba4", "v"),  # magenta
    "445": ("#008300", "P"),  # green
    "555": ("#eda100", "X"),  # yellow
    "935": ("#e34948", "*"),  # red
}
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------
def list_view_files(data_dir, pattern=DEFAULT_PATTERN):
    """Return the L1C files in ``data_dir`` sorted forward -> nadir -> aft."""
    files = glob.glob(os.path.join(data_dir, pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern!r} in {data_dir}")

    def gimbal(path):
        with h5py.File(path, "r") as h5:
            g = _attr(h5[FILE_ATTRS], "Gimbal angle average (degrees)", np.nan)
        return -g if np.isfinite(g) else 0.0

    return sorted(files, key=gimbal)


def view_tag(path):
    """Extract the view label (e.g. '646F') from an AirMSPI file name."""
    m = re.search(r"_(\d{3}[FAN])_", os.path.basename(path))
    return m.group(1) if m else os.path.basename(path)


def _attr(group, name, default=None):
    """Read a scalar/string attribute, decoding bytes and unwrapping arrays."""
    if name not in group.attrs:
        return default
    v = group.attrs[name]
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode()
    if isinstance(v, np.ndarray) and v.size == 1:
        v = v.item()
    if isinstance(v, np.generic):
        v = v.item()
    return v


def available_bands(h5):
    """Band labels (e.g. ['355', ..., '935']) present under HDFEOS/GRIDS."""
    bands = []
    for name in h5[GRID_ROOT]:
        m = re.fullmatch(r"(\d+)nm_band", name)
        if m:
            bands.append(m.group(1))
    return sorted(bands, key=int)


def is_polarized(h5, band):
    return "Q_scatter" in h5[f"{GRID_ROOT}/{band}nm_band/Data Fields"]


def channel_info(h5):
    """Map channel name (e.g. '660I') -> (center wavelength nm, F0 at 1 AU)."""
    ci = h5["Channel_Information"]
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in ci["Channel_name"][()]]
    wl = ci["Center_wavelength"][()]
    f0 = ci["Solar_irradiance_at_1_AU"][()]
    return {n: (float(w), float(e)) for n, w, e in zip(names, wl, f0)}


def grid_coords(h5):
    """UTM easting (XDim) and northing (YDim) of the pixel centers, in m."""
    return h5[f"{GRID_ROOT}/XDim"][()], h5[f"{GRID_ROOT}/YDim"][()]


def utm_zone(h5):
    proj = h5[f"{GRID_ROOT}/Ancillary/Data Fields/UTM_projection"]
    return int(_attr(proj, "utm_zone_number"))


# ---------------------------------------------------------------------------
# Location -> grid index
# ---------------------------------------------------------------------------
def utm_to_index(h5, easting, northing):
    """Nearest grid (x, y) index for a UTM easting/northing (same zone)."""
    xdim, ydim = grid_coords(h5)
    ix = int(round((easting - xdim[0]) / (xdim[1] - xdim[0])))
    iy = int(round((northing - ydim[0]) / (ydim[1] - ydim[0])))
    if not (0 <= ix < xdim.size and 0 <= iy < ydim.size):
        raise ValueError(
            f"UTM ({easting}, {northing}) is outside the grid "
            f"(E {xdim.min()}-{xdim.max()}, N {ydim.min()}-{ydim.max()})"
        )
    return ix, iy


def latlon_to_index(h5, lat, lon):
    """Nearest grid (x, y) index for a latitude/longitude (needs pyproj)."""
    try:
        from pyproj import Transformer
    except ImportError as exc:
        raise ImportError("--latlon requires pyproj; use --x/--y or --utm instead") from exc
    epsg = (32600 if lat >= 0 else 32700) + utm_zone(h5)
    tr = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    easting, northing = tr.transform(lon, lat)
    return utm_to_index(h5, easting, northing)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def _window(ix, iy, box, nx, ny):
    """Row/column slices of a box x box window centered on (ix, iy)."""
    if box < 1 or box % 2 == 0:
        raise ValueError("--box must be a positive odd integer")
    if not (0 <= ix < nx and 0 <= iy < ny):
        raise ValueError(f"Pixel (x={ix}, y={iy}) outside grid of {nx} x {ny}")
    h = box // 2
    return slice(max(iy - h, 0), min(iy + h + 1, ny)), slice(max(ix - h, 0), min(ix + h + 1, nx))


def _read(fields, name, ys, xs):
    """Read a window of a field as float64 with fill values set to NaN."""
    if name not in fields:
        return None
    a = fields[name][ys, xs].astype(np.float64)
    a[a <= FILL_VALUE] = np.nan
    return a


def _nanmean(a, valid):
    if a is None:
        return np.nan
    v = a[valid & np.isfinite(a)]
    return float(v.mean()) if v.size else np.nan


def _circmean_deg(a, valid):
    v = a[valid & np.isfinite(a)]
    if v.size == 0:
        return np.nan
    r = np.radians(v)
    return float(np.degrees(np.arctan2(np.sin(r).mean(), np.cos(r).mean())) % 360.0)


def signed_view_zenith(vza, vaa, heading, gimbal):
    """VZA signed + for forward views and - for aft views.

    View_azimuth points from the target toward the sensor, so for a forward
    view it is roughly opposite to the aircraft heading.
    """
    if not np.isfinite(vza):
        return np.nan
    c = np.cos(np.radians(vaa - heading)) if np.isfinite(vaa) and np.isfinite(heading) else np.nan
    if np.isfinite(c) and abs(c) > 0.05:
        sign = -np.sign(c)
    elif np.isfinite(gimbal) and gimbal != 0:
        sign = np.sign(gimbal)
    else:
        sign = 1.0
    return float(sign * vza)


def read_pixel(files, ix, iy, bands=None, box=1, frame="scatter", units="brf"):
    """Read multi-angle I, Q, U, DoLP and geometry for one grid pixel.

    Parameters
    ----------
    files : list of str
        L1C files, one per viewing angle.
    ix, iy : int
        Grid column (x) and row (y) index.
    bands : list of str or None
        Band labels, e.g. ['470', '660', '865']. None reads all bands.
    box : int
        Odd window size; values are averaged over valid pixels in the
        box x box window centered on (ix, iy).
    frame : {'scatter', 'meridian'}
        Reference frame for Q and U.
    units : {'brf', 'radiance'}
        'brf' applies pi * d^2 / (F0 * cos(SZA)) to I, Q and U;
        'radiance' returns W m-2 sr-1 nm-1.

    Returns
    -------
    pandas.DataFrame
        One row per (view, band).
    """
    rows = []
    for path in files:
        with h5py.File(path, "r") as h5:
            fa = h5[FILE_ATTRS]
            sun_dist = float(_attr(fa, "Sun distance", 1.0))
            heading = float(_attr(fa, "Aircraft heading (degrees)", np.nan))
            gimbal = float(_attr(fa, "Gimbal angle average (degrees)", np.nan))
            epoch_str = _attr(fa, "Epoch (UTC)", None)
            epoch = None
            if epoch_str:
                epoch = datetime.strptime(epoch_str.rstrip("Z"), "%Y-%m-%dT%H:%M:%S.%f")

            chan = channel_info(h5)
            xdim, ydim = grid_coords(h5)
            ys, xs = _window(ix, iy, box, xdim.size, ydim.size)
            cy, cx = iy - ys.start, ix - xs.start  # center inside window

            anc = h5[f"{GRID_ROOT}/Ancillary/Data Fields"]
            lat = float(anc["Latitude"][iy, ix])
            lon = float(anc["Longitude"][iy, ix])
            elev = float(anc["Elevation"][iy, ix])

            band_list = available_bands(h5) if bands is None else bands
            for band in band_list:
                grp = f"{GRID_ROOT}/{band}nm_band/Data Fields"
                if grp not in h5:
                    raise KeyError(f"Band {band} nm not found in {os.path.basename(path)}")
                f = h5[grp]
                pol = is_polarized(h5, band)

                I = _read(f, "I", ys, xs)
                imask = f["I.mask"][ys, xs] == 1
                valid_i = imask & np.isfinite(I)

                sza = _read(f, "Sun_zenith", ys, xs)
                wl, f0 = chan.get(f"{band}I", (float(band), np.nan))
                if units == "brf":
                    factor = np.pi * sun_dist**2 / (f0 * np.cos(np.radians(sza)))
                else:
                    factor = np.ones_like(I)

                rec = {
                    "view": view_tag(path),
                    "band": band,
                    "wavelength_nm": wl,
                    "polarized": pol,
                    "I": _nanmean(I * factor, valid_i),
                    "Q": np.nan,
                    "U": np.nan,
                    "DOLP": np.nan,
                    "n_valid_I": int(valid_i.sum()),
                    "n_valid_QU": 0,
                }

                if pol:
                    Q = _read(f, f"Q_{frame}", ys, xs)
                    U = _read(f, f"U_{frame}", ys, xs)
                    valid_qu = (
                        valid_i
                        & (f["Q.mask"][ys, xs] == 1)
                        & (f["U.mask"][ys, xs] == 1)
                        & np.isfinite(Q)
                        & np.isfinite(U)
                    )
                    i_qu = _nanmean(I * factor, valid_qu)
                    q = _nanmean(Q * factor, valid_qu)
                    u = _nanmean(U * factor, valid_qu)
                    rec["Q"], rec["U"] = q, u
                    # DoLP from the (box-averaged) Stokes parameters; for
                    # box=1 this equals the file's DOLP field.
                    rec["DOLP"] = float(np.hypot(q, u) / i_qu) if i_qu > 0 else np.nan
                    rec["n_valid_QU"] = int(valid_qu.sum())

                vza = _nanmean(_read(f, "View_zenith", ys, xs), valid_i)
                vaa = _circmean_deg(_read(f, "View_azimuth", ys, xs), valid_i)
                saa = _circmean_deg(_read(f, "Sun_azimuth", ys, xs), valid_i)
                t = _read(f, "Time_in_seconds_from_epoch", ys, xs)
                t_c = t[cy, cx] if t is not None else np.nan

                rec.update(
                    {
                        "Scattering_angle": _nanmean(_read(f, "Scattering_angle", ys, xs), valid_i),
                        "View_zenith": vza,
                        "VZA_signed": signed_view_zenith(vza, vaa, heading, gimbal),
                        "View_azimuth": vaa,
                        "Sun_zenith": _nanmean(sza, valid_i),
                        "Sun_azimuth": saa,
                        # Both azimuths use the photon-propagation convention.
                        "Relative_azimuth": (vaa - saa) % 360.0 if np.isfinite(vaa + saa) else np.nan,
                        "Glint_angle": _nanmean(_read(f, "Glint_angle", ys, xs), valid_i),
                        "RDQI": _nanmean(_read(f, "RDQI", ys, xs), valid_i),
                        "time_utc": (epoch + timedelta(seconds=float(t_c))).isoformat()
                        if epoch is not None and np.isfinite(t_c)
                        else "",
                        "x": ix,
                        "y": iy,
                        "utm_x": float(xdim[ix]),
                        "utm_y": float(ydim[iy]),
                        "latitude": lat,
                        "longitude": lon,
                        "elevation_m": elev,
                        "box": box,
                        "units": units,
                        "frame": frame,
                        "file": os.path.basename(path),
                    }
                )
                rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def _style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRIDLINE, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.xaxis.label.set_color(INK_SECONDARY)
    ax.yaxis.label.set_color(INK_SECONDARY)


def plot_pixel(df, out=None):
    """Plot I, Q, U, DoLP vs signed view zenith (top) and scattering angle (bottom).

    Points in the scattering-angle row are not connected: neighbouring
    scattering angles can come from very different view geometries.
    Returns the matplotlib Figure.
    """
    import matplotlib.pyplot as plt

    units = df["units"].iloc[0]
    frame = df["frame"].iloc[0]
    unit_lbl = "BRF" if units == "brf" else r"W m$^{-2}$ sr$^{-1}$ nm$^{-1}$"
    order = sorted(df["band"].unique(), key=int)
    pol_bands = [b for b in order if df.loc[df["band"] == b, "polarized"].any()]
    # Legend order: polarimetric bands first, then intensity-only bands.
    bands = pol_bands + [b for b in order if b not in pol_bands]

    panels = [
        ("I", f"I  [{unit_lbl}]", bands),
        ("Q", f"Q ({frame} plane)  [{unit_lbl}]", pol_bands),
        ("U", f"U ({frame} plane)  [{unit_lbl}]", pol_bands),
        ("DOLP", "DoLP", pol_bands),
    ]
    xaxes = [
        ("VZA_signed", "Signed view zenith angle [deg]  (+ forward, - aft)", True),
        ("Scattering_angle", "Scattering angle [deg]", False),
    ]

    fig, axes = plt.subplots(
        len(xaxes), len(panels), figsize=(17, 8), sharex="row", sharey="col",
        layout="constrained", facecolor=SURFACE,
    )
    for r, (xcol, xlabel, connect) in enumerate(xaxes):
        for c, (var, ylabel, band_set) in enumerate(panels):
            ax = axes[r, c]
            _style_axes(ax)
            if var in ("Q", "U"):
                ax.axhline(0.0, color=BASELINE, linewidth=0.9, zorder=1)
            for band in band_set:
                sub = df[df["band"] == band].dropna(subset=[xcol, var]).sort_values(xcol)
                if sub.empty:
                    continue
                color, marker = BAND_STYLE.get(band, (INK_MUTED, "o"))
                ax.plot(
                    sub[xcol], sub[var], color=color, marker=marker,
                    linestyle="-" if connect else "none", linewidth=1.5,
                    markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.0,
                    label=f"{band} nm" if band in pol_bands else f"{band} nm (I only)",
                    zorder=3,
                )
            if not band_set:
                ax.text(0.5, 0.5, "no polarimetric band selected", transform=ax.transAxes,
                        ha="center", va="center", color=INK_MUTED, fontsize=9)
            ax.set_xlabel(xlabel, fontsize=9)
            ax.set_ylabel(ylabel, fontsize=9)
            if r == 0:
                ax.set_title(var if var != "DOLP" else "DoLP", color=INK_PRIMARY, fontsize=11, loc="left")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="outside lower center", ncol=len(labels), frameon=False,
        fontsize=10, labelcolor=INK_SECONDARY,
    )

    first = df.iloc[0]
    views = ", ".join(df.drop_duplicates("view")["view"])
    fig.suptitle(
        f"AirMSPI L1C  pixel x={first['x']}, y={first['y']}  "
        f"(UTM {first['utm_x']:.0f} E, {first['utm_y']:.0f} N; "
        f"lat {first['latitude']:.5f}, lon {first['longitude']:.5f}; "
        f"elev {first['elevation_m']:.0f} m; box {first['box']}x{first['box']})\n"
        f"views: {views}",
        color=INK_PRIMARY, fontsize=11,
    )

    if out:
        fig.savefig(out, dpi=150, facecolor=SURFACE)
        print(f"Saved figure: {out}")
    return fig


def quicklook(files, ix, iy, band="660", step=8, out=None):
    """Map of the nadir (or middle) view I with the selected pixel marked.

    Shaded pixels lie inside the reference-view footprint but are missing in
    at least one other view (outside the common overlap, or masked, e.g. by
    terrain occlusion at oblique angles). Unshaded pixels have full
    multi-angle coverage. Returns the matplotlib Figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib import patheffects
    from matplotlib.colors import ListedColormap

    nadir = [p for p in files if view_tag(p).endswith("N")]
    ref = nadir[0] if nadir else files[len(files) // 2]
    grp = f"{GRID_ROOT}/{band}nm_band/Data Fields"

    all_valid = None
    for path in files:
        with h5py.File(path, "r") as h5:
            m = h5[f"{grp}/I.mask"][::step, ::step] == 1
        all_valid = m if all_valid is None else (all_valid & m)
    with h5py.File(ref, "r") as h5:
        ny, nx = h5[f"{grp}/I"].shape
        img = h5[f"{grp}/I"][::step, ::step].astype(np.float64)
    img[img <= FILL_VALUE] = np.nan
    ref_valid = np.isfinite(img)
    partial = ref_valid & ~all_valid

    # Sample k sits at full-resolution index k*step.
    h = step / 2.0
    extent = (-h, (img.shape[1] - 1) * step + h, (img.shape[0] - 1) * step + h, -h)

    fig, ax = plt.subplots(figsize=(8, 8), layout="constrained", facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    vmax = np.nanpercentile(img, 99) if ref_valid.any() else 1.0
    im = ax.imshow(img, cmap="gray", vmin=0, vmax=vmax, extent=extent, interpolation="nearest")
    ax.imshow(
        np.where(partial, 1.0, np.nan), cmap=ListedColormap([BAND_STYLE["660"][0]]),
        vmin=0, vmax=1, alpha=0.55, extent=extent, interpolation="nearest",
    )
    ax.plot(
        ix, iy, marker="+", color=BAND_STYLE["470"][0], markersize=22, markeredgewidth=2.5,
        path_effects=[patheffects.withStroke(linewidth=5, foreground=SURFACE)],
    )

    # Crop to the reference footprint (and the pixel, if outside it).
    rows, cols = np.nonzero(ref_valid)
    if rows.size:
        pad = 25 * step
        ax.set_xlim(max(min(cols.min() * step, ix) - pad, 0), min(max(cols.max() * step, ix) + pad, nx))
        ax.set_ylim(min(max(rows.max() * step, iy) + pad, ny), max(min(rows.min() * step, iy) - pad, 0))

    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.set_xlabel("x (column index)", color=INK_SECONDARY)
    ax.set_ylabel("y (row index)", color=INK_SECONDARY)
    ax.set_title(
        f"{view_tag(ref)} {band} nm I.  Shaded: missing in at least one of the "
        f"{len(files)} views.  +: pixel ({ix}, {iy})",
        fontsize=10, color=INK_PRIMARY, loc="left",
    )
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label(r"I [W m$^{-2}$ sr$^{-1}$ nm$^{-1}$]", color=INK_SECONDARY)
    cb.outline.set_edgecolor(BASELINE)
    cb.ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    if out:
        fig.savefig(out, dpi=120, facecolor=SURFACE)
        print(f"Saved quicklook: {out}")
    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Plot AirMSPI L1C multi-angle I, Q, U for one grid pixel.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    loc = p.add_argument_group("pixel location (give one)")
    loc.add_argument("--x", type=int, help="grid column index (XDim)")
    loc.add_argument("--y", type=int, help="grid row index (YDim)")
    loc.add_argument("--utm", type=float, nargs=2, metavar=("EASTING", "NORTHING"),
                     help="UTM coordinates in m (zone of the files)")
    loc.add_argument("--latlon", type=float, nargs=2, metavar=("LAT", "LON"),
                     help="latitude/longitude in degrees (needs pyproj)")

    p.add_argument("--data-dir", default=DEFAULT_DATA_DIR, help="directory with L1C files")
    p.add_argument("--pattern", default=DEFAULT_PATTERN, help="glob pattern for L1C files")
    p.add_argument("--bands", nargs="+", default=["all"],
                   help="bands in nm (e.g. 445 660), 'all', or 'pol' for the polarimetric "
                        "bands only; Q, U and DoLP exist only for polarimetric bands")
    p.add_argument("--box", type=int, default=1, help="odd window size for averaging")
    p.add_argument("--units", choices=("brf", "radiance"), default="brf")
    p.add_argument("--frame", choices=("scatter", "meridian"), default="scatter",
                   help="reference frame for Q and U")
    p.add_argument("--out", help="output figure path (default: airmspi_pixel_x<X>_y<Y>.png)")
    p.add_argument("--csv", help="also write the extracted values to this CSV file")
    p.add_argument("--quicklook", nargs="?", const="", default=None, metavar="PATH",
                   help="also save a map with the pixel marked "
                        "(default path: airmspi_quicklook_x<X>_y<Y>.png)")
    p.add_argument("--show", action="store_true", help="open interactive plot windows")
    args = p.parse_args(argv)

    n_loc = sum([args.x is not None or args.y is not None, args.utm is not None, args.latlon is not None])
    if n_loc != 1:
        p.error("give exactly one of --x/--y, --utm or --latlon")
    if (args.x is None) != (args.y is None):
        p.error("--x and --y must be given together")
    if args.box < 1 or args.box % 2 == 0:
        p.error("--box must be a positive odd integer")
    return args


def main(argv=None):
    args = parse_args(argv)

    import matplotlib

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = list_view_files(args.data_dir, args.pattern)

    with h5py.File(files[0], "r") as h5:
        if args.utm is not None:
            ix, iy = utm_to_index(h5, *args.utm)
        elif args.latlon is not None:
            ix, iy = latlon_to_index(h5, *args.latlon)
        else:
            ix, iy = args.x, args.y
        all_bands = available_bands(h5)
        pol_bands = [b for b in all_bands if is_polarized(h5, b)]

    sel = [b.lower().replace("nm", "") for b in args.bands]
    bands = all_bands if sel == ["all"] else pol_bands if sel == ["pol"] else sel
    missing = [b for b in bands if b not in all_bands]
    if missing:
        sys.exit(f"Unknown band(s) {missing}; available: {all_bands}")

    print(f"Reading {len(files)} views at pixel x={ix}, y={iy} (box {args.box}) ...")
    df = read_pixel(files, ix, iy, bands=bands, box=args.box, frame=args.frame, units=args.units)

    bad = df.loc[df["n_valid_I"] == 0, ["view", "band"]]
    if not bad.empty:
        print("Warning: no valid I for:", ", ".join(f"{v}/{b}" for v, b in bad.itertuples(index=False)))
        print("         use --quicklook to see where all views overlap")

    cols = ["view", "band", "VZA_signed", "Scattering_angle", "Sun_zenith", "Relative_azimuth",
            "I", "Q", "U", "DOLP", "n_valid_I"]
    with pd.option_context("display.width", 200, "display.float_format", "{:.5f}".format):
        print(df[cols].to_string(index=False))

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"Saved CSV: {args.csv}")

    out = args.out or f"airmspi_pixel_x{ix}_y{iy}.png"
    plot_pixel(df, out=out)

    if args.quicklook is not None:
        ql = args.quicklook or f"airmspi_quicklook_x{ix}_y{iy}.png"
        quicklook(files, ix, iy, out=ql)

    if args.show:
        plt.show()
    plt.close("all")


if __name__ == "__main__":
    main()
