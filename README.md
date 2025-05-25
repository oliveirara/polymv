# polymv

<a href="http://ascl.net/2007.009"><img src="https://img.shields.io/badge/ascl-2007.009-blue.svg?colorB=262255" alt="ascl:2007.009" /></a> [![GitHub license](https://img.shields.io/github/license/oliveirara/polymv)](https://github.com/oliveirara/polymv/blob/master/LICENSE) ![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/oliveirara/polymv)


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
- GCC or another default Unix compiler
- [Miniforge](https://github.com/conda-forge/miniforge) package manager

> 🛑 **Note:** This package was developed for Unix-based systems. Windows support is not guaranteed. \
> \
> The installation steps were designed to use `conda` due to the large amount of binaries available.

### Steps

1. Clone the repository
   ```bash
   git clone https://github.com/oliveirara/polymv.git
   cd polymv
   ```

2. Review the [Makefile](./Makefile)
   It will install the following user-level dependencies:
   - MPSolve
   - chealpix

3. Customize the [config.mk](./config.mk)
   Adjust these settings as needed:
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
For more insight into what the installer is doing, check out the [scripts](./scripts/) folder.

## 🧪 Usage

Explore the [notebooks](./notebooks/) folder for hands-on examples and demos. These notebooks will guide you through converting `alms` and analyzing cosmological data using `polymv`.

## 📝 Citation

If you use `polymv` in your research, please cite:

```bibtex
@article{2020PDU....3000608O,
       author = {{Oliveira}, Renan A. and {Pereira}, Thiago S. and {Quartin}, Miguel},
        title = "{CMB statistical isotropy confirmation at all scales using multipole vectors}",
      journal = {Physics of the Dark Universe},
     keywords = {Observational cosmology, Cosmic microwave background, Statistical isotropy, Multipole vectors, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, General Relativity and Quantum Cosmology},
         year = 2020,
        month = dec,
       volume = {30},
          eid = {100608},
        pages = {100608},
          doi = {10.1016/j.dark.2020.100608},
archivePrefix = {arXiv},
       eprint = {1812.02654},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020PDU....3000608O},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

and

```bibtex
@article{2024arXiv241108087R,
       author = {{Rodrigues}, Ricardo G. and {Pereira}, Thiago S. and {Quartin}, Miguel},
        title = "{Fr{\'e}chet Vectors as sensitive tools for blind tests of CMB anomalies}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics, General Relativity and Quantum Cosmology},
         year = 2024,
        month = nov,
          eid = {arXiv:2411.08087},
        pages = {arXiv:2411.08087},
          doi = {10.48550/arXiv.2411.08087},
archivePrefix = {arXiv},
       eprint = {2411.08087},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv241108087R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

## 💰 Funding

This work was supported by:

- Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq)
- Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES)
- Fundação Araucária (PBA-2016)

## 📄 License

Licensed under the GNU General Public License v3.
See the [LICENSE](./LICENSE) file for details.
