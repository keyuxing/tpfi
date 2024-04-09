# tpfi
[![PyPI version](https://badge.fury.io/py/tpfi.svg)](https://badge.fury.io/py/tpfi)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`tpfi` is an easy-to-use visualization tool for astronomers to identify and analyze 
stars in the Kepler, K2, and TESS missions. The main focus of this project is on 
two functions: `plot_identification` and `plot_season`. These functions create 
plots to help visualize target stars and their surrounding environments.

---

**Plot identification charts for Kepler, K2 and TESS.**

[![pPKxsjH.png](https://s1.ax1x.com/2023/08/14/pPKxsjH.png)](https://imgse.com/i/pPKxsjH)

[![pPKxcDA.png](https://s1.ax1x.com/2023/08/14/pPKxcDA.png)](https://imgse.com/i/pPKxcDA)

[![pPKx6ud.png](https://s1.ax1x.com/2023/08/14/pPKx6ud.png)](https://imgse.com/i/pPKx6ud)

The `plot_identification` function creates identification charts, which are useful 
for determining if a target is contaminated by nearby stars. In each chart, the 
right panel overlays the Gaia DR3 catalog onto the TESS Target Pixel Files (TPF) 
with the target marked by a cross symbol. The circle size represents the relative 
brightness of the stars according to Gaia G magnitude. The left panel displays the 
same sky coverage but taken from the 
[DSS2 Red survey](https://skyview.gsfc.nasa.gov/current/cgi/moreinfo.pl?survey=DSS2%20Red).

This function is revised based on 
[_tpfplotter_](https://github.com/jlillo/tpfplotter). 

---

**Plot season charts for Kepler targets.**

[![pPKxrge.png](https://s1.ax1x.com/2023/08/14/pPKxrge.png)](https://imgse.com/i/pPKxrge)

The Kepler Space Telescope was designed to observe a specific field of stars 
continuously, but it also needed to ensure its solar panels faced the Sun to power 
its operations. As Kepler orbited the Sun, it performed a 90-degree roll every three 
months (or every quarter) to keep its solar panels sunward. After four of these 
quarterly rolls, amounting to a full year and a complete 360-degree rotation, Kepler 
would return to its original orientation. Consequently, every set of quarters 
separated by four (e.g., Q1, Q5, Q9) would find Kepler with the same pointing towards 
its target field, as the telescope effectively revisited its initial orientation after 
its annual orbit.

The `plot_season` function creates a plot of the TPF for each season of the Kepler 
mission for a given target star. Since the aperture mask may change even in the same 
season, the function only plots the aperture mask for the first quarter of each
season. This function is only applicable for Kepler targets.

Season 1: Q1, Q5, Q9, Q13, Q17

Season 2: Q2, Q6, Q10, Q14

Season 3: Q3, Q7, Q11, Q15

Season 4: Q4, Q8, Q12, Q16

## Installation

You can install this package using `pip`:
```shell
pip install tpfi
```

To upgrade to the latest version, run:
```shell
pip install -U tpfi
```

## How to use

See the [example notebook](https://github.com/keyuxing/tpfi/blob/main/examples/tutorial.ipynb) for more details.

## How to cite

If you use this package in your research, please add the following sentence to the footnote 
or acknowledgement section of your paper:
> This work has made use of \texttt{tpfi} (publicly available at \url{https://github.com/keyuxing/tpfi}).

and cite our paper: 

Xing et al., 2024 [[ads link](https://ui.adsabs.harvard.edu/abs/2024ApJS..271...57X)][[iop link](https://doi.org/10.3847/1538-4365/ad2ddd)]

## Contributing

If you would like to contribute to this project, feel free to submit a pull request 
or open an issue on GitHub. Any suggestion, improvement, or bug report is welcomed.

## License

This project is licensed under the MIT License - see the 
[LICENSE](https://github.com/keyuxing/tpfi/blob/main/LICENSE) file for details.