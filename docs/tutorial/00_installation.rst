Installation
============

RNFtools is a complex package based on Python. Its installation should be done with PIP_ but there exist
also other options. Since all neccesary programs are compiled
automatically, some libraries and tools are required to be present in your system.

.. contents::
	:depth: 3


Requirements
------------

Requirements for basic installation of RNFtools are:

* Unix-like operating system (Linux, MacOS, etc.).
* `Python`_ 3.2+.


Additional requirements
^^^^^^^^^^^^^^^^^^^^^^^

RNFtools installs all required programs on fly when they are requested. However, some libraries are required for successful compilation and if a required library is missing, a problem will not occur during installation of RNFtools but during execution of RNFtools scripts.

**On Linux**

* *GCC 4.7+*
* *zlib* library.
* *pdflib* library (it is needed for GNUplot PDF outputs). `PDFlib Lite`_ can be used.
* `GIT`_

**On OSX**

* *XCode*.
* *pdflib* library (it is needed for GNUplot PDF outputs). `PDFlib Lite`_ can be used (``brew install pdflib-lite``).
* `GIT`_


Installation with root account
------------------------------

Installation using PIP from PyPI (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On most of machines, RNFtools can be installed using `PIP`_ from `PyPI`_ by 

.. code-block:: bash
	
	pip3 install rnftools


If this command does not work, check if PIP3 is installed in your system (the command may have a slightly different name, e.g., ``pip``, ``pip-3``, ``pip3.4``, ``pip-3.4``). If not, install PIP by the `official instructions`_ (or try ``easy_install3 pip``).

Upgrade to the newest version can be done also using `PIP`_.

.. code-block:: bash

	pip3 install --upgrade rnftools


Installation using Easy Install from PyPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install RNFtools also using `Easy Install`_:

.. code-block:: bash

	easy_install3 rnftools


Installation using PIP from GIT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to install RNFtools directly from `GIT repository`_, enter these commands:

.. code-block:: bash

	git clone git://github.com/karel-brinda/rnftools
	pip3 install rnftools


Installation without root account
---------------------------------

If you do not have root permissions, the simplest solution is to install some Python distribution
into your home directory. Our recommendation is `Anaconda`_. Then you can install RNFtools exactly
in the same way as described in previous section.

If you insist on your main Python installation, you have to use update few variables and then use again
the same procedure as in the previous section.


First create a directory where RNFtools will be installed.

.. code-block:: bash
	
	mkdir ~/rnftools


Then save this directory into variable ``PYTHONUSERBASE`` 

.. code-block:: bash
	
	export PYTHONUSERBASE=~/rnftools


Now you can install RNFtools. The parameter ``--user`` implies installation into the predefined directory. 

.. code-block:: bash
	
	pip3 install --user rnftools


As the last step, add these lines into your ``~/.bashrc``


.. code-block:: bash

	export PYTHONUSERBASE=~/rnftools
	export PATH=$PATH:~/rnftools/bin




.. _`official instructions`: https://pip.pypa.io/en/latest/installing.html
.. _`GIT`: https://git-scm.com/
.. _`Python`: https://www.python.org/downloads/
.. _`Anaconda`: http://continuum.io/downloads
.. _`SnakeMake`: http://bitbucket.org/johanneskoester/snakemake/
.. _`PIP`: http://pip.pypa.io/en/latest/installing.html
.. _`PyPI`: https://pypi.python.org/pypi
.. _`Easy Install`: http://pypi.python.org/pypi/setuptools
.. _`GIT repository`: http://github.com/karel-brinda/rnftools
.. _`PDFlib lite`: http://www.pdflib.com/download/free-software/pdflib-lite-7/