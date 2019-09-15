# GeospatialML101

## This repository was written and tested on a Debian Ubuntu Linux machine running PopOS


## Linux Libspatialindex binaries prerequisite installation:
```bash
geopandas.sjoin() requires rtree with libspatialindex binaries. In the project 
main directory, execute the following procedure.

https://github.com/libspatialindex/libspatialindex/wiki/1.-Getting-Started

This tutorial will help you to get started with libspatialindex using C++ on linux. 
The following code is run on Ubuntu Disco. If you are using windows the installation 
may be different. First install some prerequisites please enter the following into 
terminal. You may very well already have these installed.

sudo apt-get update
sudo apt-get install curl
sudo apt-get install g++
sudo apt-get install make

Now we download and install the library. It doesn't matter what directory you 
download this in. Please note the version number, you can check if there are 
more recent versions in the download page here:
 http://download.osgeo.org/libspatialindex/ . 

Now enter the following into your terminal (cd to project directory first!):

curl -L http://download.osgeo.org/libspatialindex/spatialindex-src-1.8.5.tar.gz | tar xz
cd spatialindex-src-1.8.5
./configure
make
sudo make install
sudo ldconfig
```

## Setup (Editable)
```bash
# Using this repository as the current working directory
python3.7 -m venv .pyenv  # Install Python 3.7 virtual environment
.pyenv/bin/pip install --upgrade pip wheel pip-tools bumpversion tox  # Install additional tools
.pyenv/bin/pip install -e .  # Install this application (in editable mode)
```

## Run the unit tests
```bash
# Using this repository as the current working directory
# Note: the dask scheduler does not print to stdout because it's running in a different process.
.pyenv/bin/tox
```

## Generate the latest dependencies (to update requirements.txt)
```bash
# Using this repository as the current working directory
.pyenv/bin/pip-compile -vvv --upgrade --dry-run setup.py
```