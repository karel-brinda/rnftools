Installation
------------

Requirements
^^^^^^^^^^^^

Requirements for basic installation of RNFtools.

* Unix-like operating system (Linux, MacOS, etc.).
* `Python`_ 3.2+.

	* If not installed yet, it can be downloaded it from https://www.python.org/downloads/.

Additional requirements
"""""""""""""""""""""""

RNFtools is installing all required programs on fly when requisted. Some needed libraries are required for successful compilation. If some of the libraries is missing, a problem will not occur during installation of RNFtools but during execution of RNFtools scripts.

**On Linux**

* *GCC 4.7+*
* *zlib* library.

**On OSX**

* *XCode*.
* *pdflib* library. It can be installed by ``brew install pdflib-lite``.
    
	* RNFtools use PNG and PDF GnuPlot terminals which cannot compiled without this library. 


Installation using PIP
^^^^^^^^^^^^^^^^^^^^^^

On most of machines, RNFtools can be installed using `PIP`_ by 

.. code-block:: bash
	
	pip install rnftools

If this command does not work, check if PIP3 is installed in your system (the command may have a slightly different name, e.g., ``pip3``, ``pip-3``, ``pip3.4``, ``pip-3.4``). If not, install PIP by instructions on https://pip.pypa.io/en/latest/installing.html and repeat the previous step.

Upgrade to the newest version can be done by `PIP`_, too.

.. code-block:: bash

	pip install --upgrade rnftools


Installation using Easy Install
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install RNFtools also using `Easy Install`_:

.. code-block:: bash

	easy_install rnftools


Installation from GIT
^^^^^^^^^^^^^^^^^^^^^

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
