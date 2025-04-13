include config.mk
SHELL := /usr/bin/bash

all: create_env install_dependencies install_gmp install_mpsolve install_cfitsio install_chealpix install_nlopt install_polymv set_env

create_env:
	@bash scripts/create_env.sh ${PKG_MANAGER} ${PYTHON_ENV} ${PYTHON_VERSION}

install_dependencies:
	@bash scripts/install_dependencies.sh ${PKG_MANAGER} ${PYTHON_ENV}

install_gmp:
	@bash scripts/install_gmp.sh ${PKG_MANAGER} ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_mpsolve:
	@bash scripts/install_mpsolve.sh ${PKG_MANAGER} ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_cfitsio:
	@bash scripts/install_cfitsio.sh ${PKG_MANAGER} ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_chealpix:
	@bash scripts/install_chealpix.sh ${PKG_MANAGER} ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_nlopt:
	@bash scripts/install_nlopt.sh ${PKG_MANAGER} ${PYTHON_ENV} ${INSTALLATION_FOLDER}

install_polymv:
	@echo "🔧 Installing polymv..."
	git clone https://github.com/oliveirara/polymv.git && \
	cd polymv && \
	conda activate ${PYTHON_ENV} && \
	python setup.py install && \
	cd .. && \
	rm -rf polymv && \
	echo "✅ polymv installed."

set_env:
	@echo "🔧 Setting environment variables..."
	export LD_LIBRARY_PATH=${INSTALLATION_FOLDER}/lib:$$LD_LIBRARY_PATH && \
	export LD_LIBRARY_PATH=${INSTALLATION_FOLDER}/lib64:$$LD_LIBRARY_PATH && \
	echo "✅ Environment variables set."
	@read -p "Do you want to add these environment variables to your .bashrc file? {y/n}: " choice; \
	if [ $$choice = y ]; then \
		echo "export LD_LIBRARY_PATH=${INSTALLATION_FOLDER}/lib:$$LD_LIBRARY_PATH" >> ~/.bashrc && \
		echo "export LD_LIBRARY_PATH=${INSTALLATION_FOLDER}/lib64:$$LD_LIBRARY_PATH" >> ~/.bashrc && \
		echo "✅ Environment variables added to .bashrc."; \
	else \
		echo "⚠️  Remember to set the environment variables manually if needed."; \
	fi