# Recovery of Microbial Genomes from Metagenomes and their functional annotation

**This workshop covers the recovery of microbial genomes from metagenomics samples using [Metawrap](https://github.com/bxlab/metaWRAP/tree/master) and their functional annotation using [METABOLIC](https://github.com/AnantharamanLab/METABOLIC)**

![workflow]{images/workflow.jpg}

## Setup

### Install miniconda3 if necessary

Details on how to install miniconda are available [here](https://docs.anaconda.com/free/miniconda/#quick-command-line-install)

Below is the quick version pasted from the above link:

```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
```

You will then need to add some extra channels (repositories to search for packages), in order to be able to create the environments:
```bash
conda config --add channels conda-forge bioconda
```

### Clone this repository to your own space on Maxwell

This will give you environment files and the scripts needed to run the analysis (which will be uploaded AFTER the workshop, to force you to write them yourself for now.)

```bash
git clone LINK
```

### Create the environments

The `metawrap.yaml` environment file is provided to create the conda environment that contains all the packages / libraries needed to conduct microbial genome assembly from metagenomic data.

The `METABOLIC.yaml` environment file is provided to create the conda enviornment that contains all the packages / libraries needed to conduct functionall annotation of MAGs.

It is recommended to install Mamba, a conda 'drop-in' replacement, before creating the environments to make it quicker:
```bash
conda install -y mamba
```

You can then create the environments using:
```bash
mamba env create -y -f metawrap.yaml
```
and
```bash
mamba env create -y -f METABOLIC.yaml
```
THIS MAY NEED TO BE SOME AS SOME SBATCH --WRAP SCRIPT TO AVOID OVERLOADING THE LOGIN NODE

You should now have the environments set up to run the analyses.

## Get the Raw Data

Run the download script to get the raw data for the analysis. Do not worry about understanding how this step works just now, as it is unlikely to apply when using your own data.

```bash
sbatch Download.sh
```

## Quality Control

Before we can proceed with any sort of analysis, we need to do some quality control on the data we have.
