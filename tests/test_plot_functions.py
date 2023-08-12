import unittest
import warnings
from pathlib import Path

import lightkurve as lk
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from parameterized import parameterized_class

from tpfi import plot_identification, plot_season, plot_sky, plot_tpf

OUTPUT_DIR = Path("test_output")
OUTPUT_DIR.mkdir(exist_ok=True)


@parameterized_class(
    [
        {"target": "TIC150428135", "search_params": {"sector": 1, "exptime": 120, "author": "SPOC"}},
        {"target": "KIC1161345", "search_params": {"quarter": 17, "exptime": 1800}},
        {"target": "EPIC211765471", "search_params": {"campaign": 5, "exptime": 1800}},
    ]
)
class TestIdentificationFunctions(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.tpf = lk.search_targetpixelfile(self.target, **self.search_params).download()

    def tearDown(self):
        plt.close("all")

    def test_plot_sky(self):
        plot_sky(self.tpf)
        plt.savefig(OUTPUT_DIR / f"test_sky_{self.tpf.mission}.png", bbox_inches="tight", dpi=300)

    def test_plot_tpf(self):
        plot_tpf(self.tpf)
        plt.savefig(OUTPUT_DIR / f"test_tpf_{self.tpf.mission}.png", bbox_inches="tight", dpi=300)

    def test_plot_tpf_with_cbar(self):
        _, ax = plt.subplots(figsize=(6, 4))
        divider = make_axes_locatable(ax)
        ax_cb = divider.append_axes("right", size="10%", pad=0.3)

        plot_tpf(self.tpf, ax=ax, ax_cb=ax_cb)
        plt.savefig(OUTPUT_DIR / f"test_tpf_with_cbar_{self.tpf.mission}.png", bbox_inches="tight", dpi=300)

    def test_plot_identification(self):
        plot_identification(self.tpf)
        plt.savefig(OUTPUT_DIR / f"test_identification_{self.tpf.mission}.png", bbox_inches="tight", dpi=300)

    def test_plot_identification_with_params(self):
        plot_identification(
            self.tpf,
            mag_limit=20,
            frame=0,
            cmap="Blues",
            c_star="k",
            c_mask="grey",
            show_label=False,
            show_ticklabels=False,
            verbose=True,
        )
        plt.savefig(
            OUTPUT_DIR / f"test_identification_with_params_{self.tpf.mission}.png", bbox_inches="tight", dpi=300
        )


class TestSeasonFunction(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_plot_season(self):
        plot_season("KIC2991403")
        plt.savefig(OUTPUT_DIR / "test_season.png", bbox_inches="tight", dpi=300)

    def test_plot_season_with_params(self):
        plot_season("KIC2991403", mag_limit=20, cmap="Blues", c_star="k", c_mask="grey", show_label=False, verbose=True)
        plt.savefig(OUTPUT_DIR / "test_season_with_params.png", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    unittest.main()
