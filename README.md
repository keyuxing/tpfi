# tpfi
[![PyPI version](https://badge.fury.io/py/tpfi.svg)](https://badge.fury.io/py/tpfi)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`tpfi` is an easy-to-use visualization tool for astronomers to identify and analyze 
stars in the Kepler, K2, and TESS missions. The main focus of this project is on 
two functions: `plot_identification` and `plot_season`. These functions create 
plots to help visualize target stars and their surrounding environments.

---

**Plot identification charts for Kepler, K2 and TESS.**

![alt text](https://raw.githubusercontent.com/keyuxing/tpfi/main/examples/kepler.png)

![alt text](https://raw.githubusercontent.com/keyuxing/tpfi/main/examples/k2.png)

![alt text](https://raw.githubusercontent.com/keyuxing/tpfi/main/examples/tess.png)

The `plot_identification` function creates identification charts, which are useful 
for determining if a target is contaminated by nearby stars. In each chart, the 
right panel displays the target (marked by a cross symbol) and nearby stars. The 
circle size indicates their relative brightness. The left panel displays the same 
sky coverage but from the [DSS2 Red survey](https://skyview.gsfc.nasa.gov/current/cgi/moreinfo.pl?survey=DSS2%20Red).

This function is revised based on 
[_tpfplotter_](https://github.com/jlillo/tpfplotter). 

---

**Plot season charts for Kepler targets.**

![alt text](https://raw.githubusercontent.com/keyuxing/tpfi/main/examples/season.png)

The `plot_season` function creates a plot of the TPF for each season of the Kepler 
mission for a given target star. Note that this function is only applicable for 
Kepler targets.

## Installation

You can install this package using `pip`:
```shell
pip install tpfi
```

## How to use

See the [example notebook](https://github.com/keyuxing/tpfi/blob/main/examples/tutorial.ipynb) for more details.

## Contributing

If you would like to contribute to this project, feel free to submit a pull request 
or open an issue on GitHub. Any suggestion, improvement, or bug report is welcomed.

## License

This project is licensed under the MIT License - see the 
[LICENSE](https://github.com/keyuxing/tpfi/blob/main/LICENSE) file for details.