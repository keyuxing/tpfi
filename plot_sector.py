from argparse import ArgumentParser
from pathlib import Path

import lightkurve as lk
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylatex import Command, Document, Section
from pylatex.package import Package
from pylatex.utils import NoEscape
from pypdf import PdfMerger

from tpfi import plot_identification


def cli():
    parser = ArgumentParser()
    parser.add_argument("-s", "--sector", type=int, default=0)

    return str(parser.parse_args().sector)


def generate_cover(pdf_path, sector_name):
    doc = Document(document_options=["a4paper"])
    doc.preamble.append(Package("authblk"))

    doc.preamble.append(
        Command(
            "title",
            NoEscape(r"Identification charts of compact pulsator candidates from TESS Sector~" + sector_name[1:]),
        )
    )
    doc.preamble.append(Command("author", options=["1", "2"], arguments="Keyu Xing"))
    doc.preamble.append(Command("author", options=["1", "2"], arguments="Weikai Zong"))
    doc.preamble.append(Command("author", options="3", arguments=NoEscape(r"St\'ephane Charpinet")))
    doc.preamble.append(
        Command(
            "affil",
            options="1",
            arguments=NoEscape(r"Department of Astronomy, Beijing Normal University, Beijing~100875, P.~R.~China;"),
        )
    )
    doc.preamble.append(
        Command(
            "affil",
            options="2",
            arguments=NoEscape(
                r"Institute for Frontiers in Astronomy and Astrophysics, Beijing Normal University, "
                r"Beijing~102206, P.~R.~China;"
            ),
        )
    )
    doc.preamble.append(
        Command(
            "affil",
            options="3",
            arguments=NoEscape(
                r"IRAP, Universit\'e de Toulouse, CNRS, UPS, CNES, 14 avenue Edouard Belin, "
                r"F-31400, Toulouse, France"
            ),
        )
    )
    doc.preamble.append(Command("date", ""))
    doc.append(NoEscape(r"\maketitle"))

    # Add stuff to the document
    with doc.create(Section("Description")):
        doc.append(
            NoEscape(
                r"These identification charts are useful for determining if a compact target is contaminated "
                r"by nearby stars. In each chart, the right panel overlays the Gaia DR3 catalog onto the TESS "
                r"Target Pixel Files (TPF) with the target marked by a cross symbol. The circle size represents "
                r"the relative brightness of the stars according to Gaia G magnitude. The left panel displays "
                r"the same sky coverage but taken from the DSS2 Red survey"
                r"\footnote{https://skyview.gsfc.nasa.gov/current/cgi/moreinfo.pl?survey=DSS2\%20Red}. "
                r"To locate a specific target, use the hotkey ``Ctrl+F'' and enter its TIC number. "
                r"These charts are generated using *tpfi*\footnote{https://github.com/keyuxing/tpfi}."
            )
        )

    doc.generate_pdf(pdf_path)


def create_figure_axes():
    fig = plt.figure(figsize=(8.27, 11.69))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace=0.2)
    axes = [plt.subplot(gs[i]) for i in range(3)]
    plt.subplots_adjust(left=0, right=0.95, top=0.95, bottom=0.05)
    return fig, axes


SECTOR = cli()

fig, axes = create_figure_axes()
sector_name = f"S{str(SECTOR).zfill(3)}"

cover_name = sector_name + "_cover"
generate_cover(Path.cwd() / cover_name, sector_name)

plot_name = sector_name + "_plot.pdf"
with PdfPages(Path.cwd() / plot_name) as pdf:
    tpf_path_list = list(Path.cwd().glob("tpfs/*.fits"))

    for i in range((len(tpf_path_list) + 2) // 3 * 3):
        fig_index = i % 3
        if fig_index == 0:
            fig, axes = create_figure_axes()

        try:
            try:
                tpf = lk.read(tpf_path_list[i])
            except lk.LightkurveError:
                tpf_filename = str(tpf_path_list[i])
                target_label = f'TIC{tpf_filename.split("-")[2].lstrip("0")}'
                target_sector = int(tpf_filename.split("-")[1][2:])
                tpf = lk.search_targetpixelfile(
                    target_label, sector=target_sector, exptime=120, author="SPOC"
                ).download()

            print("-----------------------")
            print(f"Num {i + 1}, TIC {tpf.targetid}")

            plot_identification(tpf, axes[fig_index])

        except IndexError:
            axes[fig_index].axis("off")

        # Save
        if fig_index == 2:
            pdf.savefig()
            plt.close()

cover_name += ".pdf"
pdf_name = sector_name + ".pdf"
pdf_merger = PdfMerger()
pdf_merger.append(str(Path.cwd() / cover_name))
pdf_merger.append(str(Path.cwd() / plot_name))
pdf_merger.write(str(Path.cwd() / pdf_name))

Path.unlink(Path.cwd() / cover_name)
Path.unlink(Path.cwd() / plot_name)
