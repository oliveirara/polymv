# polymv

<a href="http://ascl.net/2007.009"><img src="https://img.shields.io/badge/ascl-2007.009-blue.svg?colorB=262255" alt="ascl:2007.009" /></a> [![GitHub license](https://img.shields.io/github/license/oliveirara/polymv)](https://github.com/oliveirara/polymv/blob/master/LICENSE) ![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/oliveirara/polymv)

Welcome to the **polymv** documentation!

`polymv` is a Python/C package for converting multipolar coefficients (`alms`) into **Multipole Vectors (MVs)** and **Fréchet Vectors (FVs)** for a given multipole. It's especially useful in cosmological studies, including the analysis of Planck 2015 and 2018 temperature maps.

🔭 You can explore the MVs and FVs derived from Planck data here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3866410.svg)](https://doi.org/10.5281/zenodo.3866410)

## Features

* Converts `alms` to Multipole Vectors (MVs) and Fréchet Vectors (FVs)
* Designed for cosmological data analysis (e.g. Planck temperature maps)
* Includes ready-to-use notebooks and utilities

## Getting Started

To get started with **polymv**, please refer to the following sections:

* [Installation](#installation)
- [Usage](./usage.md)
* [Examples](./test-example.md)

## Installation

### Prerequisites

Make sure the following are installed on your system:

* Python ≥ 3.10
* GCC or another default Unix compiler
* [Miniforge](https://github.com/conda-forge/miniforge) package manager

> 🛑 **Note:** This package was developed for Unix-based systems. Windows support is not guaranteed.

> The installation steps were designed to use `conda` due to the large amount of binaries available.

### Steps

1. Clone the repository
   ```bash
   git clone https://github.com/oliveirara/polymv.git
   cd polymv
   ```

2. Review the [Makefile](https://github.com/oliveirara/polymv/blob/main/Makefile). It will install the following user-level dependencies:
    - MPSolve
    - chealpix

3. Customize the [config.mk](https://github.com/oliveirara/polymv/blob/main/config.mk). Adjust these settings as needed:
    - Installation directory (defaults to user home directory)
    - Package manager (`mamba` by default, can be changed to `conda`)
    - Environment name
    - Python version

4. Run the installer
   ```bash
   make install
   ```

    Or changing the parameters online (no need to rewrite the `config.mk` file). For the example, I will change the environment name and Python version:
   ```bash
   make install PYTHON_ENV="<environment-name>" PYTHON_VERSION="<custom-python-version>"
   ```

That’s it! 🎉

For more insight into what the installer is doing, check out the [scripts](https://github.com/oliveirara/polymv/tree/main/scripts) folder.

## License

This project is licensed under the GNU General Public License v3. See the [LICENSE](LICENSE.md) file for details.

## Acknowledgments

This work was supported by:

- Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq)
- Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES)
- Fundação Araucária (PBA-2016)
