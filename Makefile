include config.mk
SHELL := /usr/bin/bash

install: create_env install_dependencies install_gmp install_mpsolve install_cfitsio install_chealpix install_nlopt install_polymv

create_env:
	@bash scripts/create_env.sh ${PKG_MANAGER} ${PYTHON_ENV} ${PYTHON_VERSION}

install_dependencies:
	@bash scripts/install_dependencies.sh ${PKG_MANAGER} ${PYTHON_ENV}

install_gmp:
	@bash scripts/install_gmp.sh ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_mpsolve:
	@bash scripts/install_mpsolve.sh ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_cfitsio:
	@bash scripts/install_cfitsio.sh ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_chealpix:
	@bash scripts/install_chealpix.sh ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_nlopt:
	@bash scripts/install_nlopt.sh ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_polymv:
	@bash scripts/install_polymv.sh ${PYTHON_ENV}
