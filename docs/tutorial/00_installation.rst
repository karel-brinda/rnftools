Installation
------------

Requirements
^^^^^^^^^^^^

Requirements for basic installation of RNFtools are:

* Unix-like operating system (Linux, MacOS, etc.).
* `Python`_ 3.2+.

Additional requirements
"""""""""""""""""""""""

RNFtools installs all required programs on fly when they are requested. However, some libraries are required for successful compilation and if a required library is missing, a problem will not occur during installation of RNFtools but during execution of RNFtools scripts.

**On Linux**

* *GCC 4.7+*
* *zlib* library.

**On OSX**

* *XCode*.
* *pdflib* library. It can be installed by ``brew install pdflib-lite``.


Installation using PIP
^^^^^^^^^^^^^^^^^^^^^^

On most of machines, RNFtools can be installed using `PIP`_ by 

.. code-block:: bash
	
	pip install rnftools

If this command does not work, check if PIP3 is installed in your system (the command may have a slightly different name, e.g., ``pip3``, ``pip-3``, ``pip3.4``, ``pip-3.4``). If not, install PIP by the `official instructions`_.

Upgrade to the newest version can be done also using `PIP`_.

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


.. _`official instructions`: https://pip.pypa.io/en/latest/installing.html
.. _Python: https://www.python.org/downloads/
.. _Anaconda: http://continuum.io/downloads
.. _SnakeMake: http://bitbucket.org/johanneskoester/snakemake/
.. _PIP: http://pip.pypa.io/en/latest/installing.html
.. _`Easy Install`: http://pypi.python.org/pypi/setuptools
.. _GIT repository: http://github.com/karel-brinda/rnftools
