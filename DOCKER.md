Quickstart using Docker
-----------------------

First install docker.  If running on Mac or Windows, you can use `docker-machine ip` to get the docker host vm's IP address.  The following code assumes you are using the `default` docker-machine image.  If otherwise, change accordingly.

```
git clone https://github.com/karel-brinda/rnftools.git
cd rnftools
docker-compose build
docker-compose up -d
open "http://`docker-machine ip default`:28888"
```

This will open a browser to an iPython Notebook complete with all of the rnftools examples.
