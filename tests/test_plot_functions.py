import unittest
import warnings
from pathlib import Path

import lightkurve as lk
import matplotlib.pyplot as plt

from tpfi import plot_identification, plot_season

OUTPUT_DIR = Path("test_output")
OUTPUT_DIR.mkdir(exist_ok=True)


class TestPlotFunctions(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

    def test_plot_identification_tess(self):
        tpf = lk.search_targetpixelfile("TIC150428135", sector=1, exptime=120, author="SPOC").download()
        plot_identification(tpf)
        plt.savefig(OUTPUT_DIR / "test_tess_identification.png", bbox_inches="tight")

    def test_plot_identification_kepler(self):
        tpf = lk.search_targetpixelfile("KIC1161345", quarter=17, exptime=1800).download()
        plot_identification(tpf)
        plt.savefig(OUTPUT_DIR / "test_kepler_identification.png", bbox_inches="tight")

    def test_plot_identification_k2(self):
        tpf = lk.search_targetpixelfile("EPIC211765471", campaign=5, exptime=1800).download()
        plot_identification(tpf)
        plt.savefig(OUTPUT_DIR / "test_k2_identification.png", bbox_inches="tight")

    def test_plot_season(self):
        plot_season("KIC10139564")
        plt.savefig(OUTPUT_DIR / "test_season.png", bbox_inches="tight")


if __name__ == "__main__":
    unittest.main()
