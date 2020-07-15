# dyngen in a Docker container

[dynverse/dyngen](https://hub.docker.com/r/dynverse/dyngen) contains all necessary packages to run dyngen from start to finish.

## Running the container
To run the container, you can use the following command.

```sh
docker run --rm -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true -v `pwd`:/home/rstudio/workdir dynverse/dyngen
```

<!-- fedora users currently need to run:
docker run --rm --ulimit="nofile=4096" -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true -v `pwd`:/home/rstudio/workdir dynverse/dyngen
-->

Keep this window open, and open up a browser and go to [127.0.0.1:8787](127.0.0.1:8787). Open up the file `example.R` for a small example on how to run a dyngen simulation.

The command can be dissected as follows.

```sh
docker run \

  # remove container after use
  --rm \
  
  # specify which port rstudio server uses
  -p 127.0.0.1:8787:8787 \
  
  # disable authentication because I'm lazy
  -e DISABLE_AUTH=true \
  
  # mount the current working directory to the rstudio home folder
  # so you will see it right away when rstudio starts
  -v `pwd`:/home/rstudio/ \
  
  # specify which container to run
  dynverse/dyngen
```

## Building the container

To rebuild this docker container from scratch, run the following command.

```sh
docker build -t dynverse/dyngen --build-arg GITHUB_PAT=$GITHUB_PAT .
```

GITHUB_PAT should be an environment variable corresponding to the Personal Access Token created by following [this tutorial](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token).


