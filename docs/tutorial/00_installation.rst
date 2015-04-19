Installation
------------

RNFtools can work on any Unix-like computer with pre-installed Python 3.2+ (if Python is not installed yet, download it from https://www.python.org/downloads/).

On most of machines, RNFtools can be installed using PIP by 

.. code-block:: bash
	
	pip install rnftools

If this command does not work, check if PIP3 is installed in your system (the command may have a slightly different name, e.g., ``pip3``, ``pip-3``, ``pip3.4``, ``pip-3.4``). If not, install PIP by instructions on https://pip.pypa.io/en/latest/installing.html and repeat the previous step.

You can install RNFtools also using Easy Install:

.. code-block:: bash

	easy_install rnftools

If you want to install RNFtools directly from GIT, use these steps:

.. code-block:: bash

	git clone git://github.com/karel-brinda/rnftools
	cd rnftools
	./install.sh

Upgrade to the newest version can be done by PIP, too.

.. code-block:: bash

	pip install --upgrade rnftools
