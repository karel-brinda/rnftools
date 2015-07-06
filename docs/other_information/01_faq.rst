.. _faq:

FAQ
===

Frequently asked questions.

.. contents::
   :depth: 3


A script which used to work does not work any more
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Try to upgrade RNFtools to the latest version

	.. code-block:: bash
	
		pip install --upgrade rnftools
	
Sometimes developers change the way how their software is compiled. Such change must be reflected in the SMBL library which is employed by RNFtools for required software installation. Upgrade of RNFtools causes upgrade of `SMBL`_ to the latest version as well.

If the problem still appears, please send us an e-mail to karel.brinda@gmail.com.


LAVEnder  
^^^^^^^^


How to upgrade a used program to the last version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to upgrade some of used programs to the latest version support by `SMBL`_, remove the program from ``~/.smbl/bin``, then it will be downloaded and installed again. When you remove the whole directory, all programs will be re-installed.


How to use a specific version of a program
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Replace the program in ``~/.smbl/bin`` by your desired version.

.. _`SMBL`: http://github.com/karel-brinda/smbl/
