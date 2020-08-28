# Experiments for StartLink

Georgia Institute of Technology, Atlanta, Georgia, USA

Reference: PAPER LINK


## Overview

This repository contains the data and source code needed to reproduce all results found in the StartLink paper.

## Program Versions

StartLink requires external software. To reproduce the results, it is recommended that you have the exact versions 
listed here, and they should be available through your environment's `$PATH`:

- ClustalO: `v1.2.4`
- DIAMOND: `0.9.24`

StartLink is a `python` based program. To get all packages, it is recommended that the user creates a `conda` environment from the file 
in `install/startlink.yml` through the following command:

    conda env create -f install/startlink.yml --name startlink

This can then be activated via 

    conda activate startlink

See `info/reproduce.[html|pdf]` for more information.

## StartLink Gene Predictions

Due to the nature of StartLink as a comparative approach, experiments can take several hours to several days to complete (depending on your setup). In light of this, we also provide all intermediary results, which include all gene predictions done by StartLink(+). These can be found in the `results` directory.


##  Reproducing Results

We provide a document detailing how to reproduce all results. This can be found at `info/reproduce.[html|pdf]`


## Folder structure

The following directory structure should be maintained while running experiments.

    .
    ├── bin                                   # Executables constructed from python/bash drivers (via install.sh)
    ├── bin_external                          # External tools
    ├── config                                # Configuration files, e.g. StartLink parameters
    ├── config.sh                             # Load bash variables for paths to different directories
    ├── db                                    # Location of constructed databases
    ├── install                               # Conda environment file for easy installation
    ├── lists                                 # Lists of genomes (main input method to scripts)
    ├── info                                  # Information about reproducing results
    ├── metadata                              # Non-genomic data, including taxonomy information
    ├── data                                  # Data Location: where all raw data will be stored during runs
    │   ├── GCFID 1                           # ID of genome 1
    │   │   ├── ncbi.gff                      # RefSeq annotation
    │   │   ├── sequence.fasta                # Genomic sequence file
    │   ├── GCFID 2                           # ID of genome 2
    │   │   ├── ncbi.gff                      # RefSeq annotation
    │   │   ├── sequence.fasta                # Genomic sequence file
    │   │   ...
    ├── code                                  # Source code
    │   ├── python                            # Python code
    │   │   ├── driver                        # Drivers that can be executed
    │   │   ├── lib                           # Library files
    │   ├── bash                              # Bash scripts
    │   │   ├── driver                        # Drivers that can be executed
    │   │   ├── lib                           # Library files
    ├── runs                                  # Data Location: where all raw data will be stored during runs
    │   ├── GCFID 1                           # ID of genome 1
    │   │   ├── startlink                     # StartLink runs
    │   │   ├── gms2                          # GMS2 runs
    │   ├── GCFID 2                           # ID of genome 2
    │   │   ├── startlink                     # StartLink runs
    │   │   ├── gms2                          # GMS2 runs
    │   │   ...
