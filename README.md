# Recovery of Microbial Genomes from Metagenomes and their functional annotation

**This workshop covers the recovery of microbial genomes from metagenomics samples using [Metawrap](https://github.com/bxlab/metaWRAP/tree/master) and their functional annotation using [METABOLIC](https://github.com/AnantharamanLab/METABOLIC)**

![workflow]{images/workflow.png}

## Setup

### Install miniconda3 if necessary

Details on how to install miniconda are available [here](https://docs.anaconda.com/free/miniconda/#quick-command-line-install)

Below is the quick version pasted from the above link. Make sure to navigate to your home directory first using `cd`.

Install miniconda:
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

Navigate to somewhere sensible
```bash
cd /uoa/scratch/users/your-username
```

```bash
git clone https://github.com/biomoff/UoA_Metagenomics_Workshop
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
sbatch --nodes=1 --cpus-per-task=1 --mem-per-cpu=4G --wrap "mamba env create -y -f metawrap.yaml > /dev/null"
```
and
```bash
sbatch --nodes=1 --cpus-per-task=1 --mem-per-cpu=4G --wrap "mamba env create -y -f METABOLIC.yaml > /dev/null"
```

*Note: we might usually just use `mamba env create -y -f file.yaml` directly, but here we are submitting it to the cluster to avoid overloading the login node, particularly if all trying to do it at the same time. We are also discarding the text that is usually shown on the terminal when create environments interactively*

**You should now have the environments set up to run the analyses.**

## Get the Raw Data

Run the download script to get the raw data for the analysis. Do not worry about understanding how this step works just now, as it is unlikely to apply when using your own data.

```bash
sbatch Scripts/Download.sh
```

You should end up with a directory called `RAW_READS` containing some metagenomics sequence files:
```
A_1.fastq
A_2.fastq
B_1.fastq
B_2.fastq
```

## Quality Control

Before we can proceed with any sort of analysis, we need to do some quality control on the data we have.
