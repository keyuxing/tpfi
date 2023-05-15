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
    tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile], verbose: bool
) -> Union[Table, None]:
    """
    Query the objects in the area of TPF from Gaia Catalog (Now Gaia DR3).

    Parameters
    ----------
    tpf: `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`
    verbose: bool
        Whether to print the Gaia Source ID of the target

    Returns
    -------
    r: `astropy.table.Table` or None
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
    r = j.get_results()

    if not (r["dist"] < cross_threshold).any():
        if verbose:
            print("Target not found in Gaia DR3")
        return None
    else:
        if verbose:
            print(f"Target Gaia Source DR3 ID: {r[0]['source_id']}")
        return r


def plot_sky(ax_sky: Axes, tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile], show_label: bool, verbose: bool):
    """
    Plot the sky image corresponding to the target pixel file.

    Parameters
    ----------
    ax_sky: `matplotlib.axes.Axes`
        A matplotlib axes object to plot the sky image into
    tpf: `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`
    show_label: bool
        Whether to show the label of the target
    verbose: bool
        Whether to show the progress of the query
    """

    # Use pixel scale for query size
    pixel_scale = 21 if tpf.meta["TELESCOP"] == "TESS" else 4

    x_pixel, y_pixel = tpf.shape[1:][1], tpf.shape[1:][0]
    x_length = x_pixel * pixel_scale
    y_length = y_pixel * pixel_scale

    theta, reverse = calculate_theta(tpf.wcs.wcs)

    # Sky
    if verbose:
        print("Querying Sky Image...")
    sky_img = query_sky_img(tpf, theta, x_length, y_length, reverse)

    ax_sky.imshow(sky_img, origin="lower")
    ax_sky.set_xticks([])
    ax_sky.set_yticks([])
    ax_sky.set_xticklabels([])
    ax_sky.set_yticklabels([])
    ax_sky.set_xlim(0, sky_img.shape[1])
    ax_sky.set_ylim(0, sky_img.shape[0])
    ax_sky.invert_xaxis()

    if show_label:
        ax_sky.add_artist(AnchoredText(tpf.meta["OBJECT"], frameon=False, loc="upper left", prop=dict(size=13)))

    # Add orientation arrows
    add_orientation(ax=ax_sky, theta=theta, pad=0.15, color="k", reverse=reverse)

    # Add scale bar
    add_scalebar(ax=ax_sky, pad=0.15, length=1 / x_pixel, scale='{}"'.format(pixel_scale))


def plot_tpf(
    ax_tpf: Axes,
    tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile],
    r: Table,
    cmap: str,
    c_star: str,
    c_mask: str,
    show_ticklabels: bool,
    mag_limit: float,
    ax_cb: Axes = None,
):
    """
    Plot the identification charts.

    Parameters
    ----------
    ax_tpf: `matplotlib.axes.Axes`
        A matplotlib axes object to plot the tpf into
    tpf: `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`
    r: `astropy.table.Table`
        The table of the Gaia objects in the area of TPF
    cmap: str
        The colormap to use
    c_star: str
        The color of the stars
    c_mask: str
        The color of the pipeline mask
    show_ticklabels: bool
        Whether to show the tick labels
    mag_limit: float
        The magnitude limit (Gaia G mag) to plot the stars in the TPF
    ax_cb: `matplotlib.axes.Axes`
        A matplotlib axes object to plot the color bar into
    """

    # TPF plot
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        median_flux = np.nanmedian(tpf.flux.value, axis=0)

    try:
        division = int(np.log10(np.nanmax(median_flux)))
    except ValueError:
        division = 0

    image = median_flux / 10**division

    im, norm = imshow_norm(image, ax_tpf, stretch=SqrtStretch(), origin="lower", cmap=cmap, zorder=0)
    x_pixel, y_pixel = tpf.shape[1:][1], tpf.shape[1:][0]
    ax_tpf.set_xlim(-0.5, x_pixel - 0.5)
    ax_tpf.set_ylim(-0.5, y_pixel - 0.5)

    if show_ticklabels:
        ax_tpf.set_xticks(np.arange(0, x_pixel, 1))
        ax_tpf.set_yticks(np.arange(0, y_pixel, 1))
        ax_tpf.set_xticklabels(np.arange(1, x_pixel + 1, 1))
        ax_tpf.set_yticklabels(np.arange(1, y_pixel + 1, 1))
    else:
        ax_tpf.set_xticks([])
        ax_tpf.set_yticks([])

    ax_tpf.yaxis.set_ticks_position("right")
    ax_tpf.invert_xaxis()

    if r is None:
        at = AnchoredText("No Gaia DR3 Data", frameon=False, loc="upper left", prop=dict(size=13))
        ax_tpf.add_artist(at)
    else:
        target_gaia_id = r[0]["source_id"]
        target_gaia_mag = r[0]["phot_g_mean_mag"]

        r.sort("phot_g_mean_mag")
        this = np.nonzero(r["source_id"] == target_gaia_id)[0][0]
        magnitude_limit = max(target_gaia_mag + 3, mag_limit)
        r = r[r["phot_g_mean_mag"] < magnitude_limit][: max(this + 50, 300)]

        qr = QTable([r["ra"].filled(), r["dec"].filled(), r["pmra"].filled(0), r["pmdec"].filled(0)])
        coords_gaia = SkyCoord(
            qr["ra"], qr["dec"], pm_ra_cosdec=qr["pmra"], pm_dec=qr["pmdec"], frame="icrs", obstime=REF_EPOCH
        )
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ErfaWarning)
            coords_obs = coords_gaia.apply_space_motion(new_obstime=tpf.time[0])

        x, y = tpf.wcs.world_to_pixel(coords_obs)
        gaia_mags = np.asarray(r["phot_g_mean_mag"])

        size_k = 1.2 * np.piecewise(
            target_gaia_mag,
            [target_gaia_mag < 12, 12 <= target_gaia_mag < mag_limit, target_gaia_mag > mag_limit],
            [70, lambda mag: 190 - mag * 10, 10],
        )
        if tpf.meta["TELESCOP"] != "TESS":
            size_k = size_k * 5
        size = size_k / 1.5 ** (gaia_mags - target_gaia_mag)

        ax_tpf.scatter(x, y, s=size, c=c_star, alpha=0.5, edgecolor=None, zorder=11)
        ax_tpf.scatter(x[this], y[this], marker="x", c="white", s=size_k / 2.5, zorder=12)

    # Pipeline aperture
    aperture = tpf.pipeline_mask
    aperture_masks = [(i, j) for i in range(aperture.shape[0]) for j in range(aperture.shape[1]) if aperture[i, j]]
    for i, j in aperture_masks:
        xy = (j - 0.5, i - 0.5)
        if aperture[i, j]:
            ax_tpf.add_patch(patches.Rectangle(xy, 1, 1, color=c_mask, fill=True, alpha=0.4))
            ax_tpf.add_patch(patches.Rectangle(xy, 1, 1, color=c_mask, fill=False, alpha=0.6, lw=1.5, zorder=9))
        else:
            ax_tpf.add_patch(patches.Rectangle(xy, 1, 1, color="gray", fill=False, alpha=0.2, lw=0.5, zorder=8))

    # Add orientation arrows
    theta, reverse = calculate_theta(tpf.wcs.wcs)
    add_orientation(ax=ax_tpf, theta=theta, pad=0.15, color="k", reverse=reverse)

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
    mag_limit: float = None,
    cmap: str = "viridis",
    c_star: str = "red",
    c_mask: str = "tomato",
    show_label: bool = True,
    show_ticklabels: bool = True,
    verbose: bool = False,
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
    mag_limit : float, optional
        The magnitude limit (Gaia G mag) to plot the stars in the TPF. If not provided, the default value will be used.
        (18 for TESS, 19.5 for Kepler/K2)
    cmap : str, optional
        The colormap to use for the TPF. Default is 'viridis'.
    c_star: str, optional
        The color of the stars in the TPF. Default is 'red'.
    c_mask: str, optional
        The color of the pipeline mask in the TPF. Default is 'tomato'.
    show_label: bool, optional
        Whether to show the label of the target in the sky image. Default is True.
    show_ticklabels : bool, optional
        Whether to show the tick labels in the TPF. Default is True.
    verbose : bool, optional
        Whether to print out progress messages. Default is False.

    Returns
    -------
    ax : `matplotlib.axes.Axes`
        The matplotlib axes object.
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4))

    if mag_limit is None:
        mag_limit = 18 if tpf.meta["TELESCOP"] == "TESS" else 19.5

    divider = make_axes_locatable(ax)
    ax_tpf = divider.append_axes("right", size="100%", pad=0.1)
    ax_cb = divider.append_axes("right", size="8%", pad=0.35)

    plot_sky(ax, tpf, show_label, verbose)

    r = query_nearby_gaia_objects(tpf, verbose=verbose)
    plot_tpf(ax_tpf, tpf, r, cmap, c_star, c_mask, show_ticklabels, mag_limit, ax_cb)


def plot_season(
    label: str,
    ax: Axes = None,
    mag_limit: float = 19.5,
    cmap: str = "gray_r",
    c_star: str = "red",
    c_mask: str = "tomato",
    show_label: bool = True,
    verbose: bool = False,
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
    mag_limit : float, optional
        The magnitude limit (Gaia G mag) to plot the stars in the TPF. Default is 19.5.
    c_star: str, optional
        The color of the stars in the TPFs. Default is 'red'.
    c_mask: str, optional
        The color of the pipeline mask in the TPFs. Default is 'tomato'.
    show_label: bool, optional
        Whether to show the label of the target in the sky image. Default is True.
    cmap : str, optional
        The colormap to use for the TPF. Default is 'gray_r'.
    verbose : bool, optional
        Whether to print out progress messages. Default is False.

    Returns
    -------
    ax : `matplotlib.axes.Axes`
        The matplotlib axes object.
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
            if quarter % 4 == season:
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

    plot_sky(ax, tpf_list[max_index], show_label, verbose)

    r = query_nearby_gaia_objects(tpf_list[max_index], verbose=verbose)
    for i, tpf in enumerate(tpf_list):
        if tpf is not None:
            plot_tpf(ax_list[i], tpf_list[i], r, cmap, c_star, c_mask, False, mag_limit)
            at = AnchoredText(f"Season {i + 1}", frameon=False, loc="upper left", prop=dict(size=13), zorder=100)
            ax_list[i].add_artist(at)
