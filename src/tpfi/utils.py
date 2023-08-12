from typing import Tuple, Union

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
from astropy.wcs import Wcsprm
from astroquery.hips2fits import hips2fits
from lightkurve import KeplerTargetPixelFile, TessTargetPixelFile
from matplotlib.axes import Axes


def calculate_theta(w: Wcsprm) -> Tuple[float, bool]:
    """
    Calculates the rotation angle of TPF according to the wcs in the FITS Header.

    Parameters
    ----------
    w: `~astropy.wcs.Wcsprm`

    Returns
    -------
    theta: float
        Rotation angle of TPF [degree]
    reverse: bool
        Whether the direction of the TPF is reversed
    """

    pc = w.pc
    cdelt = w.cdelt
    cd = cdelt * pc

    det = np.linalg.det(cd)
    sgn = np.sign(det)

    theta = -np.arctan2(sgn * cd[0, 1], sgn * cd[0, 0])
    cdelt1 = sgn * np.sqrt(cd[0, 0] ** 2 + cd[0, 1] ** 2)

    if cdelt1 < 0:
        return theta + np.pi, True
    else:
        return theta, False


def add_orientation(ax: Axes, theta: float, pad: float, color: str, reverse: bool):
    """
    Plot the orientation arrows.

    Parameters
    ----------
    ax: `matplotlib.axes.Axes`
        A matplotlib axes object to plot into
    theta: float
        Rotation angle of TPF [degree]
    pad: float
        The padding between the arrow base and the edge of the axes
    color: str
        The color of the orientation arrows
    reverse: bool
        Whether the direction of the TPF is reversed
    """

    def get_arrow_loc():
        return pad * np.cos(theta) * 0.45 * ratio, pad * np.sin(theta) * 0.45

    def get_text_loc(x, y):
        return 1 - (pad * ratio + 1.7 * x), pad + 1.7 * y

    ratio = ax.get_data_ratio()
    x1, y1 = get_arrow_loc()
    theta += -np.pi / 2 if reverse else np.pi / 2
    x2, y2 = get_arrow_loc()

    ax.arrow((1 - pad * ratio), pad, -x1, y1, color=color, head_width=0.01, transform=ax.transAxes, zorder=100)
    ax.text(
        *get_text_loc(x1, y1),
        s="E",
        color=color,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10,
        zorder=100
    )

    ax.arrow((1 - pad * ratio), pad, -x2, y2, color=color, head_width=0.01, transform=ax.transAxes, zorder=100)
    ax.text(
        *get_text_loc(x2, y2),
        s="N",
        color=color,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10,
        zorder=100
    )


def add_scalebar(ax: Axes, pad: float, length: float, scale: str):
    """
    Plot the scale bar.

    Parameters
    ----------
    ax: `matplotlib.axes.Axes`
        A matplotlib axes object to plot into
    pad: float
        The padding between the left endpoint of the scale bar and the edge of the axes
    length: float
        The length of the scale bar
    scale: str
        Scale of the scale bar
    """

    ratio = ax.get_data_ratio()
    x_min = pad * ratio / 2
    x_max = x_min + length
    y_min = pad * 0.95
    y_max = pad * 1.05
    ax.hlines(y=pad, xmin=x_min, xmax=x_max, colors="k", ls="-", lw=1.5, transform=ax.transAxes)
    ax.vlines(x=x_min, ymin=y_min, ymax=y_max, colors="k", ls="-", lw=1.5, transform=ax.transAxes)
    ax.vlines(x=x_max, ymin=y_min, ymax=y_max, colors="k", ls="-", lw=1.5, transform=ax.transAxes)
    ax.text(x_min + length / 2, pad * 1.1, scale, horizontalalignment="center", fontsize=10, transform=ax.transAxes)


def query_sky_img(
    tpf: Union[TessTargetPixelFile, KeplerTargetPixelFile],
    theta: float,
    x_length: float,
    y_length: float,
    reverse: bool,
    verbose: bool,
) -> np.ndarray:
    """
    Query the image of the area of TPF from DSS2 Red Survey.

    Parameters
    ----------
    tpf: `lightkurve.TessTargetPixelFile` or `lightkurve.KeplerTargetPixelFile`
        Target pixel files read by `lightkurve`
    theta: float
        Rotation angle of TPF [degree]
    x_length: float
        The x length of the TPF [arcsecond]
    y_length: float
        The y length of the TPF [arcsecond]
    reverse: bool
        Whether the direction of the TPF is reversed
    verbose: bool
        Whether to show the progress of querying sky image.

    Returns
    -------
    2-D array
    """

    def query_sky_data(hips):
        return hips2fits.query(
            hips=hips,
            width=n_pixel,
            height=n_pixel,
            projection="TAN",
            ra=center_ra * u.deg,
            dec=center_dec * u.deg,
            fov=radius * u.arcsec,
            format="png",
            rotation_angle=Angle(-theta * u.rad),
            cmap="Greys",
        )

    n_pixel = 200 if tpf.meta["TELESCOP"] == "TESS" else 50

    center_ra, center_dec = tpf.wcs.all_pix2world([(tpf.shape[1:][1] + 1) / 2], [(tpf.shape[1:][0] + 1) / 2], 1)

    if reverse:
        theta = np.pi - theta

    radius = 1.5 * max(x_length, y_length)

    try:
        if verbose:
            print("Querying Sky Image from DSS2 Red...")
        sky_data = query_sky_data("CDS/P/DSS2/red")
    except Exception:
        if verbose:
            print("Querying Sky Image from DSS2 Red failed. Retry with DSS2 NIR...")
        sky_data = query_sky_data("CDS/P/DSS2/NIR")

    if reverse:
        sky_data = np.flip(sky_data, axis=0)
    else:
        sky_data = np.flip(sky_data, axis=(0, 1))

    x_pixel_sky = n_pixel * x_length / radius
    x_start = int((n_pixel - x_pixel_sky) / 2)
    x_stop = n_pixel - x_start

    y_pixel_sky = n_pixel * y_length / radius
    y_start = int((n_pixel - y_pixel_sky) / 2)
    y_stop = n_pixel - y_start

    return sky_data[y_start:y_stop, x_start:x_stop]
