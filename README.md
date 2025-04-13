# polymv

<a href="http://ascl.net/2007.009"><img src="https://img.shields.io/badge/ascl-2007.009-blue.svg?colorB=262255" alt="ascl:2007.009" /></a> [![GitHub license](https://img.shields.io/github/license/oliveirara/polymv)](https://github.com/oliveirara/polymv/blob/master/LICENSE)

`polymv` is a Python/C package for converting multipolar coefficients (`alms`) into **Multipole Vectors (MVs)** and **Fréchet Vectors (FVs)** for a given multipole. It's especially useful in cosmological studies, including the analysis of Planck 2015 and 2018 temperature maps.

🔭 You can explore the MVs and FVs derived from Planck data here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3866410.svg)](https://doi.org/10.5281/zenodo.3866410)

---

## ✨ Features

- Converts `alms` to Multipole Vectors (MVs) and Fréchet Vectors (FVs)
- Designed for cosmological data analysis (e.g. Planck temperature maps)
- Includes ready-to-use notebooks and utilities

## ⚙️ Installation

### Prerequisites

Make sure the following are installed on your system:

- Python ≥ 3.10
- GCC and GFortran compilers
- [Conda or Mamba](https://github.com/conda-forge/miniforge) package manager

> 🛑 **Note:** This package was developed for Unix-based systems. Windows support is not guaranteed.

### Steps

1. Clone the repository
   ```bash
   git clone https://github.com/oliveirara/polymv.git
   cd polymv
   ```

2. Review the [Makefile](./Makefile)
   It will install the following system-level dependencies:
   - GMP
   - MPSolve
   - CFITSIO
   - chealpix
   - NLOPT

3. Customize the [config.mk](./config.mk)  
   Adjust these settings as needed:
   - Installation directory (defaults to your home directory)
   - Package manager (`mamba` by default, can be changed to `conda`)
   - Environment name
   - Python version

4. Run the installer
   ```bash
   make install
   ```

That’s it! 🎉  
For more insight into what the installer is doing, check out the [scripts](./scripts/) folder.

## 🧪 Usage

Explore the [notebooks](./notebooks/) folder for hands-on examples and demos. These notebooks will guide you through converting `alms` and analyzing cosmological data using `polymv`.

## 📝 Citation

If you use `polymv` in your research, please cite:

> R. A. Oliveira, T. S. Pereira, and M. Quartin,  
> **CMB statistical isotropy confirmation at all scales using multipole vectors**,  
> [Phys. Dark Univ. 30 (2020) 100608](https://doi.org/10.1016/j.dark.2020.100608)  
> ([arXiv:1812.02654](https://arxiv.org/abs/1812.02654))

## 💰 Funding

This work was supported by:

- Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq)  
- Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES)  
- Fundação Araucária (PBA-2016)

## 📄 License

Licensed under the GNU General Public License v3.  
See the [LICENSE](./LICENSE) file for details.
