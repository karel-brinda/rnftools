RNFtools
========

.. image:: https://travis-ci.org/karel-brinda/rnftools.svg?branch=master
	:target: https://travis-ci.org/karel-brinda/rnftools

.. image:: https://readthedocs.org/projects/rnftools/badge/?version=latest
	:target: http://rnftools.rtfd.org


Read Naming Format is a generic format for assigning
read names with encoded information about original positions. RNFtools is an associated
software package which can:

* simulate RNF-compliant reads using a wide class of integrated read simulators (Art, DwgSim, Mason, WgSim, etc.);
* evaluate mappers using RNF reads;
* convert non-RNF simulated reads to RNF (e.g., from SAM format);
* transform genomic coordinates of RNF reads between different coordinate systems (using chain LiftOver format).

Quickstart using Docker
-----------------------

First install docker.  If running on Mac or Windows, you can use `docker-machine ip` to get the docker host vm's IP address.  The following code assumes you are using the `default` docker-machine image.  If otherwise, change accordingly.

```
git clone https://github.com/karel-brinda/rnftools.git
cd rnftools
docker-compose build
docker-compose up -d

```

This will open a browser to an iPython Notebook complete with all of the rnftools examples.

Links
-----

**Web of the project:** http://karel-brinda.github.io/rnftools/

**RNF specification:** http://karel-brinda.github.io/rnf-spec/

**Documentation:** http://rnftools.rtfd.org

**Installation:** http://rnftools.readthedocs.org/en/latest/tutorial/00_installation.html

**Examples of usage:** http://github.com/karel-brinda/rnftools/tree/master/examples/tutorial

**Publication:** http://dx.doi.org/10.1093/bioinformatics/btv524 (preprint: http://arxiv.org/abs/1504.00556)
