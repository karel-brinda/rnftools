.PHONY: \
	all clean install \
	test test_repo test_repo_coverage test_parallel test_package \
	pylint flake8 yapf coverage \
	inc pypi sha256 \
	docs readme wpypi wconda \
	deppip depconda \
	submodules \
	help

PIP=/usr/bin/env pip
PYTHON=/usr/bin/env python3

ROOT_DIR = $(shell pwd)

###############
# BASIC RULES #
###############

all: rnftools ## Compile RNFtools

help: ## Print help message
	    @echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	$(PYTHON) setup.py clean --all
	rm -fr _index_test/ _test_*
	$(MAKE) -C tests clean
	($(MAKE) -C docs clean || true) > /dev/null 2> /dev/null

install: ## Install RNFtools using PIP
	$(PIP) uninstall -y rnftools || true
	$(PIP) install .

###########
# TESTING #
###########

test: test_repo

coverage: ## Run test coverage analysis
coverage: test_repo_coverage

test_repo: ## Run unit tests & integration from the repo dir
test_repo:
	$(MAKE) -C tests clean
	$(MAKE) -C tests

test_repo_coverage: ## Compute test coverage
	$(MAKE) -C tests clean
	# replace /usr/bin/env/python3 by coverage
	PATH=$$(pwd)/bin/python3_coverage_wrapper:$$PATH $(MAKE) -C tests

test_parallel: ## Run tests in parallel
	$(MAKE) -C tests clean
	$(MAKE) -C tests parallel

pylint: ## Run PyLint
	$(PYTHON) -m pylint rnftools

flake8: ## Run Flake8
	flake8

yapf: ## Run YAPF (inline replacement)
	yapf -i --recursive rnftools setup.py tests


#############
# RELEASING #
#############

inc: ## Increment version
	./increment_version.py

pypi: ## Upload RNFtools to PyPI
	$(MAKE) clean
	$(PYTHON) setup.py sdist bdist_wheel upload

sha256: ## Compute sha256 for the PyPI package
sha256:
	s=$$(curl https://pypi.python.org/pypi/rnftools  2>/dev/null| perl -pe 's/#/\n/g' | grep -o 'https.*\.tar\.gz' | xargs curl -L 2>/dev/null | shasum -a 256 | awk '{print $$1;}'); echo $$s; echo $$s | pbcopy


#######################
# DOCUMENTATION & WEB #
#######################

docs: ## Build and open Sphinx documentation
	$(MAKE) -C docs html
	open docs/.build/html/index.html || true

readme: ## Convert README to HTML
	rst2html README.rst > README.html

wconda: ## Open RNFtools Bioconda webpage
	open https://bioconda.github.io/recipes/rnftools/README.html

wpypi: ## Open RNFtools PyPI webpage
	open https://pypi.python.org/pypi/rnftools


########################
# INSTALL DEPENDENCIES #
########################

depconda: ## Install dependencies using Conda
	cat requirements.txt | perl -pe 's/==.*//g' | xargs conda install

deppip: ## Install dependencies using PIP
	cat requirements.txt | xargs $(PIP) install

