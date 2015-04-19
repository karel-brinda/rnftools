Installation
------------

Requirements
""""""""""""

* Unix-like operating system (Linux, MacOS, etc.).
* `Python`_ 3.2+. If not installed yet, download it from https://www.python.org/downloads/. This minimum version is required by `SnakeMake`_ on which RNFtools are based.

Installation using PIP
""""""""""""""""""""""

On most of machines, RNFtools can be installed using `PIP`_ by 

.. code-block:: bash
	
	pip install rnftools

If this command does not work, check if PIP3 is installed in your system (the command may have a slightly different name, e.g., ``pip3``, ``pip-3``, ``pip3.4``, ``pip-3.4``). If not, install PIP by instructions on https://pip.pypa.io/en/latest/installing.html and repeat the previous step.

Upgrade to the newest version can be done by `PIP`_, too.

.. code-block:: bash

	pip install --upgrade rnftools


Installation using Easy Install
"""""""""""""""""""""""""""""""

You can install RNFtools also using `Easy Install`_:

.. code-block:: bash

	easy_install rnftools


Installation from GIT
"""""""""""""""""""""

If you want to install RNFtools directly from `GIT repository`_, enter these commands:

.. code-block:: bash

	git clone git://github.com/karel-brinda/rnftools
	cd rnftools
	./install.sh

.. _Python: http://python.org
.. _Anaconda: http://continuum.io/downloads
.. _SnakeMake: http://bitbucket.org/johanneskoester/snakemake/
.. _PIP: http://pip.pypa.io/en/latest/installing.html
.. _`Easy Install`: http://pypi.python.org/pypi/setuptools
.. _GIT repository: http://github.com/karel-brinda/rnftools
