# DES Quick Release

This is a codebase for generating the DES quick release catalogs (i.e., Y2Q1, Y3Q2, etc.). The basic idea is to make catalog-level coadds based on the single epoch images, thus circumventing the process of building image level coadds.

## Installation

There is not much support provided for installation. However, everything is in python using standard packages, so cloning is a good start:
```
git clone https://github.com/kadrlica/desqr.git
```

## Execution

This code base is essentially a pipeline for processing single epoch catalogs produced by DESDM and converting them into a catalog-level coadd. The pipeline is executed in appropriately named stages:
```
./run_01.0_download.py config.yaml 
./run_02.0_pixelize.py config.yaml 

...                   

./run_07.0_release.py  config.yaml 
```

Each pipeline stage has its own argument parser which can be accessed through the conventional `-h` or `--help` command line argument. Each pipeline stage takes the same yaml configuration file (`config.yaml`) as input.
