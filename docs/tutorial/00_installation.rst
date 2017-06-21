Installation
============

RNFtools is distributed as a Python-based package, which is distributed through
BioConda_ (a bioinformatics channel for the Conda_ package manager) and PIP_.
Since BioConda_ does not require a root account and installs also all the
RNFtools dependencies, it is the recommended way of installation.

.. contents::
	:depth: 3


Requirements
------------

Requirements for basic installation of RNFtools are:

* Unix-like operating system (Linux, OS X, etc.).
* `Python`_ 3.3+.

When RNFtools is installed using BioConda_, all the additional dependencies are
installed automatically. If not, installation of the following programs is up
to user.

* `Art`_
* `CuReSim`_
* `DWGsim`_
* `Mason`_
* `WGsim`_
* `SamTools`_



Installation using Bioconda (recommended)
-----------------------------------------

The easiest and safest approach of RNFtools installation is to create a
separate Bioconda_ environment.

.. code-block:: bash

        conda create -n rnftools rnftools

Once the environment is installed, you can activate it by

.. code-block:: bash

        source activate rnftools

and deactivate by

.. code-block:: bash

        source deactivate


Alternatively, RNFtools can be installed directly in the default Conda_
environment.  However, this approach might not work in certain situations due
to possible collisions with the dependencies of your already installed
packages.

.. code-block:: bash

        conda install rnftools


Installation using PIP from PyPI
--------------------------------

RNFtools can be installed using `PIP`_ from `PyPI`_ by

.. code-block:: bash

	pip install rnftools


If this command does not work, check if pip is installed in your system (the
command may have a slightly different name, e.g., ``pip``, ``pip-3``,
``pip.4``, ``pip-3.4``). If not, install PIP by the `official instructions`_
(or try ``easy_install3 pip``).

Upgrade to the newest version can be done also using `PIP`_.

.. code-block:: bash

	pip install --upgrade rnftools


Installation using PIP from GIT
-------------------------------

To install RNFtools directly from `GIT repository <http://github.com/karel-brinda/rnftools>`_, run

.. code-block:: bash

	git clone git://github.com/karel-brinda/rnftools
	pip install rnftools

or

.. code-block:: bash

	pip install git+http://github.com/karel-brinda/rnftools


Installation using PIP without a root account
---------------------------------------------

First, we need to create a directory where RNFtools will be installed.

.. code-block:: bash

	mkdir ~/rnftools


Then we have to add its path into the variable ``PYTHONUSERBASE``

.. code-block:: bash

	export PYTHONUSERBASE=~/rnftools


Now we can finally install RNFtools. The parameter ``--user`` implies installation
into the predefined directory.

.. code-block:: bash

	pip install --user rnftools


As the last step, we need to add the following lines to ``~/.bashrc``

.. code-block:: bash

	export PYTHONUSERBASE=~/rnftools
	export PATH=$PATH:~/rnftools/bin


.. _`official instructions`: https://pip.pypa.io/en/latest/installing.html
.. _`Python`: https://www.python.org
.. _`Conda`: https://conda.io/
.. _`Bioconda`: https://bioconda.github.io/
.. _`SnakeMake`: https://snakemake.readthedocs.io
.. _`SamTools`: http://www.htslib.org/
.. _`PIP`: http://pip.pypa.io
.. _`PyPI`: https://pypi.python.org/pypi

.. _`Art`: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
.. _`CuReSim`: http://www.pegase-biosciences.com/curesim-a-customized-read-simulator/
.. _`DWGsim`: https://github.com/nh13/DWGSIM
.. _`Mason`: http://publications.imp.fu-berlin.de/962/
.. _`WGsim`: https://github.com/lh3/wgsim
