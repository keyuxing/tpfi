from argparse import ArgumentParser
from pathlib import Path

import requests


def cli():
    parser = ArgumentParser()
    parser.add_argument("-s", "--sector", type=int, default=0)

    return str(parser.parse_args().sector)


SECTOR = cli()
TPF_LINKS = requests.get(
    f"https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{SECTOR}_tp.sh"
).text.split("\n")[1:]
DOWNLOAD_SCRIPT = open(Path.cwd() / "download.sh", "w")
DOWNLOAD_SCRIPT.write("#!/bin/sh\n")

# Generating download script for targets of WG8
print("Generating download script...")
with open(Path.cwd() / "wg8-targets-all-latest.txt", "r") as target_list:
    for line in target_list.readlines()[10:]:
        observed_sector = line[127:362].replace(" ", "").replace("-", "").replace("U", "S").split("S")[1:]

        if SECTOR.zfill(2) in observed_sector:
            target_id = line[11:22].strip()
            for link in TPF_LINKS:
                if target_id.zfill(16) in link:
                    DOWNLOAD_SCRIPT.write("{}\n".format(link))
                    break

DOWNLOAD_SCRIPT.close()
