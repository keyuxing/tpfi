import logging
import warnings
from typing import Union

import lightkurve as lk
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Table
from astropy.time import Time
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import imshow_norm
from astroquery.gaia import Gaia
from erfa import ErfaWarning
from lightkurve import KeplerTargetPixelFile, LightkurveWarning, TessTargetPixelFile
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .utils import add_orientation, add_scalebar, calculate_theta, query_sky_img

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
Gaia.ROW_LIMIT = -1

REF_EPOCH = Time("J2016")

logging.getLogger("astroquery").setLevel(logging.WARNING)


def query_nearby_gaia_objects(
    tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile], verbose: bool = False
) -> Union[Table, None]:
    """
    Query the objects in the area of TPF from Gaia Catalog (Now Gaia DR3).

    Parameters
    ----------
    tpf: `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`
    verbose: bool, optional
        Whether to print the Gaia Source ID of the target

    Returns
    -------
    gaia_table: `astropy.table.Table` or None
        The table of the Gaia objects in the area of TPF
    """

    pixel_scale = 21 if tpf.meta["TELESCOP"] == "TESS" else 4
    tpf_radius = np.max(tpf.shape[1:]) * pixel_scale

    # Get the parameters of the Gaia sources
    radius = u.Quantity(tpf_radius * 1.5, u.arcsec)
    ra = tpf.ra * u.deg
    dec = tpf.dec * u.deg
    pm_ra = tpf.meta["PMRA"]
    pm_dec = tpf.meta["PMDEC"]
    if pm_ra and pm_dec:
        cross_threshold = 3 / 3600
        pm_ra *= u.mas / u.yr
        pm_dec *= u.mas / u.yr
        coords_j2000 = SkyCoord(ra, dec, pm_ra_cosdec=pm_ra, pm_dec=pm_dec, frame="icrs", obstime=Time("J2000"))
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ErfaWarning)
            coords = coords_j2000.apply_space_motion(new_obstime=REF_EPOCH)
    else:
        cross_threshold = 5 / 3600
        coords = SkyCoord(ra, dec, frame="icrs")

    j = Gaia.cone_search_async(coords, radius, columns=["source_id", "phot_g_mean_mag", "ra", "dec", "pmra", "pmdec"])
    gaia_table = j.get_results()

    if not (gaia_table["dist"] < cross_threshold).any():
        if verbose:
            print("Target not found in Gaia DR3")
        return None
    else:
        if verbose:
            print(f"Target Gaia Source DR3 ID: {gaia_table[0]['source_id']}")
        return gaia_table


def plot_sky(
    tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile],
    ax: Axes = None,
    show_label: bool = True,
    verbose: bool = False,
):
    """
    Plot the sky image corresponding to the target pixel file.

    Parameters
    ----------
    tpf: `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`.
    ax: `matplotlib.axes.Axes`
        A matplotlib axes object to plot the sky image into. If not provided, a new axes object will be created.
    show_label: bool, optional
        Whether to show the label of the target. Default is True.
    verbose: bool, optional
        Whether to show the progress querying sky image. Default is False.
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4))

    # Use pixel scale for query size
    pixel_scale = 21 if tpf.meta["TELESCOP"] == "TESS" else 4

    x_pixel, y_pixel = tpf.shape[1:][1], tpf.shape[1:][0]
    x_length = x_pixel * pixel_scale
    y_length = y_pixel * pixel_scale

    theta, reverse = calculate_theta(tpf.wcs.wcs)

    # Sky
    sky_img = query_sky_img(tpf, theta, x_length, y_length, reverse, verbose)

    ax.imshow(sky_img, origin="lower")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlim(0, sky_img.shape[1])
    ax.set_ylim(0, sky_img.shape[0])
    ax.invert_xaxis()

    if show_label:
        ax.add_artist(AnchoredText(tpf.meta["OBJECT"], frameon=False, loc="upper left", prop=dict(size=13)))

    # Add orientation arrows
    add_orientation(ax=ax, theta=theta, pad=0.15, color="k", reverse=reverse)

    # Add scale bar
    add_scalebar(ax=ax, pad=0.15, length=1 / x_pixel, scale='{}"'.format(pixel_scale))


def plot_tpf(
    tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile],
    ax: Axes = None,
    mag_limit: float = None,
    cmap: str = "viridis",
    c_star: str = "red",
    c_mask: str = "tomato",
    show_ticklabels: bool = True,
    ax_cb: Axes = None,
    gaia_table: Table = None,
    verbose: bool = False,
):
    """
    Plot the identification charts.

    Parameters
    ----------
    tpf: `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`.
    ax: `matplotlib.axes.Axes`
        A matplotlib axes object to plot the tpf into.
    mag_limit: float, optional
        The magnitude limit (Gaia G mag) to plot the stars in the TPF. If not provided, the default value will be used
        (18 for TESS, 19.5 for Kepler/K2).
    cmap: str, optional
        The colormap to use. Default is 'viridis'.
    c_star: str, optional
        The color of the stars. Default is 'red'.
    c_mask: str, optional
        The color of the pipeline mask. Default is 'tomato'.
    show_ticklabels: bool, optional
        Whether to show the tick labels. Default is True.
    ax_cb: `matplotlib.axes.Axes`, optional
        A matplotlib axes object to plot the color bar into. If not provided, no color bar will be plotted.
    gaia_table: `astropy.table.Table`, optional
        The table of the Gaia objects in the area of TPF. If not provided, the Gaia objects will be queried.
    verbose: bool, optional
        Whether to show the progress of querying Gaia objects. Default is False.
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4))

    if gaia_table is None:
        gaia_table = query_nearby_gaia_objects(tpf, verbose=verbose)

    if mag_limit is None:
        mag_limit = 18 if tpf.meta["TELESCOP"] == "TESS" else 19.5

    # TPF plot
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        median_flux = np.nanmedian(tpf.flux.value, axis=0)

    try:
        division = int(np.log10(np.nanmax(median_flux)))
    except ValueError:
        division = 0

    image = median_flux / 10**division

    im, norm = imshow_norm(image, ax, stretch=SqrtStretch(), origin="lower", cmap=cmap, zorder=0)
    x_pixel, y_pixel = tpf.shape[1:][1], tpf.shape[1:][0]
    ax.set_xlim(-0.5, x_pixel - 0.5)
    ax.set_ylim(-0.5, y_pixel - 0.5)

    if show_ticklabels:
        ax.set_xticks(np.arange(0, x_pixel, 1))
        ax.set_yticks(np.arange(0, y_pixel, 1))
        ax.set_xticklabels(np.arange(1, x_pixel + 1, 1))
        ax.set_yticklabels(np.arange(1, y_pixel + 1, 1))
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    ax.yaxis.set_ticks_position("right")
    ax.invert_xaxis()

    if gaia_table is None:
        at = AnchoredText("No Gaia DR3 Data", frameon=False, loc="upper left", prop=dict(size=13))
        ax.add_artist(at)
    else:
        target_gaia_id = gaia_table[0]["source_id"]
        target_gaia_mag = gaia_table[0]["phot_g_mean_mag"]

        gaia_table.sort("phot_g_mean_mag")
        this = np.nonzero(gaia_table["source_id"] == target_gaia_id)[0][0]
        magnitude_limit = max(target_gaia_mag + 3, mag_limit)
        gaia_table = gaia_table[gaia_table["phot_g_mean_mag"] < magnitude_limit][: max(this + 50, 300)]

        qr = QTable(
            [
                gaia_table["ra"].filled(),
                gaia_table["dec"].filled(),
                gaia_table["pmra"].filled(0),
                gaia_table["pmdec"].filled(0),
            ]
        )
        coords_gaia = SkyCoord(
            qr["ra"], qr["dec"], pm_ra_cosdec=qr["pmra"], pm_dec=qr["pmdec"], frame="icrs", obstime=REF_EPOCH
        )
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ErfaWarning)
            coords_obs = coords_gaia.apply_space_motion(new_obstime=tpf.time[0])

        x, y = tpf.wcs.world_to_pixel(coords_obs)
        gaia_mags = np.asarray(gaia_table["phot_g_mean_mag"])

        size_k = 1.2 * np.piecewise(
            target_gaia_mag,
            [target_gaia_mag < 12, 12 <= target_gaia_mag < mag_limit, target_gaia_mag > mag_limit],
            [70, lambda mag: 190 - mag * 10, 10],
        )
        if tpf.meta["TELESCOP"] != "TESS":
            size_k = size_k * 5
        size = size_k / 1.5 ** (gaia_mags - target_gaia_mag)

        ax.scatter(x, y, s=size, c=c_star, alpha=0.5, edgecolor=None, zorder=11)
        ax.scatter(x[this], y[this], marker="x", c="white", s=size_k / 2.5, zorder=12)

    # Pipeline aperture
    aperture = tpf.pipeline_mask
    aperture_masks = [(i, j) for i in range(aperture.shape[0]) for j in range(aperture.shape[1]) if aperture[i, j]]
    for i, j in aperture_masks:
        xy = (j - 0.5, i - 0.5)
        if aperture[i, j]:
            ax.add_patch(patches.Rectangle(xy, 1, 1, color=c_mask, fill=True, alpha=0.4))
            ax.add_patch(patches.Rectangle(xy, 1, 1, color=c_mask, fill=False, alpha=0.6, lw=1.5, zorder=9))
        else:
            ax.add_patch(patches.Rectangle(xy, 1, 1, color="gray", fill=False, alpha=0.2, lw=0.5, zorder=8))

    # Add orientation arrows
    theta, reverse = calculate_theta(tpf.wcs.wcs)
    add_orientation(ax=ax, theta=theta, pad=0.15, color="k", reverse=reverse)

    if ax_cb is not None:
        # Add color bar
        cb = Colorbar(ax=ax_cb, mappable=im, orientation="vertical", ticklocation="right")

        max_diff = int((np.nanmax(image) - np.nanmin(image)) / 0.01)
        n_ticks = min(4, max_diff)
        cbar_ticks = np.linspace(np.nanmin(image), np.nanmax(image), n_ticks, endpoint=True)

        cb.set_ticks(cbar_ticks)
        cb.set_ticklabels(np.round(cbar_ticks, 2))
        exponent = r"$\times$ $10^{}$".format(division)
        cb.set_label(r"Flux {} (e$^-$)".format(exponent), fontsize=13)


def plot_identification(
    tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile],
    ax: Axes = None,
    **kwargs,
):
    """
    Plot the identification chart for a given target pixel file (TPF).

    The function creates a combined identification chart, including the sky image and
    the TPF with stars' positions, magnitudes, and pipeline aperture.

    Parameters
    ----------
    tpf : `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`
    ax : `matplotlib.axes.Axes`, optional
        A matplotlib axes object to plot the identification chart into. If not provided, a new figure will be created.
    kwargs : dict, optional
        Other keyword arguments to be passed to `plot_sky` and `plot_tpf`.
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4))

    divider = make_axes_locatable(ax)
    ax_tpf = divider.append_axes("right", size="100%", pad=0.1)
    ax_cb = divider.append_axes("right", size="8%", pad=0.35)

    plot_sky(tpf, ax, **kwargs)
    plot_tpf(tpf, ax_tpf, ax_cb=ax_cb, **kwargs)


def plot_season(
    label: str,
    ax: Axes = None,
    verbose: bool = False,
    **kwargs,
):
    """
    Plot identification charts for different seasons of an astronomical object using data from the Kepler mission.

    The function searches for the target pixel files (TPFs) from the Kepler mission for different seasons
    and creates identification charts for each season, showing the sky image and TPFs.

    Parameters
    ----------
    label : str
        The label of the astronomical object to be searched for in the Kepler mission.
    ax : `matplotlib.axes.Axes`, optional
        A matplotlib axes object to plot the identification charts into. If not provided, a new figure will be created.
    verbose : bool, optional
        Whether to print out progress messages. Default is False.
    kwargs : dict, optional
        Other keyword arguments to be passed to `plot_sky` and `plot_tpf`.
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(18, 3))

    search_result = lk.search_targetpixelfile(label, exptime=1800, author="Kepler")
    if not len(search_result):
        raise ValueError("No TPFs found for this target.")

    quarter_array = np.full(4, -1, dtype=int)
    for season in range(4):
        for mission in search_result.mission:
            quarter = int(mission[-2:])
            if quarter != 0 and (quarter - 1) % 4 == season:
                quarter_array[season] = quarter
                break

    n_seasons = np.nonzero(quarter_array >= 0)[0].size
    if verbose:
        print(f"Found {n_seasons} seasons")

    tpf_list = []
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=LightkurveWarning)
        i = 0
        for quarter in quarter_array:
            if quarter < 0:
                tpf_list.append(None)
                continue
            if verbose:
                print(f"Downloading TPF ({i+1}/{n_seasons}) (Quarter {quarter})")
            tpf = lk.search_targetpixelfile(label, exptime=1800, quarter=quarter, author="Kepler").download()
            tpf_list.append(tpf)
            i += 1

    length_array = np.zeros(len(tpf_list))
    for i, tpf in enumerate(tpf_list):
        if tpf is not None:
            length_array[i] = tpf.shape[2]
    max_index = length_array.argmax()

    divider = make_axes_locatable(ax)
    percent_array = np.full(4, np.nan)
    for i, tpf in enumerate(tpf_list):
        if tpf is not None:
            percent_array[i] = int(tpf.shape[2] / tpf_list[max_index].shape[2] * 100)
    percent_array = np.nan_to_num(percent_array, nan=np.nanmedian(percent_array))

    ax_list = []
    for percent in percent_array:
        ax_list.append(divider.append_axes("right", size=f"{percent}%", pad=0.1))

    plot_sky(tpf_list[max_index], ax, **kwargs)

    gaia_table = query_nearby_gaia_objects(tpf_list[max_index], verbose=verbose)
    for i, tpf in enumerate(tpf_list):
        if tpf is not None:
            plot_tpf(tpf_list[i], ax_list[i], cmap="gray_r", gaia_table=gaia_table, show_ticklabels=False, **kwargs)
            at = AnchoredText(f"Season {i + 1}", frameon=False, loc="upper left", prop=dict(size=13), zorder=100)
            ax_list[i].add_artist(at)
