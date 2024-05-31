# Recovery of Microbial Genomes from Metagenomes and their functional annotation

**This workshop covers the recovery of microbial genomes from metagenomics samples using [Metawrap](https://github.com/bxlab/metaWRAP/tree/master) and their functional annotation using [METABOLIC](https://github.com/AnantharamanLab/METABOLIC)**

![workflow](images/workflow.png)

---

# Setup

## Install miniconda3 if necessary

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

If this is the first time you have used conda on Maxwell, you will need to **restart your session by logging out and back in.** 

### Adding channels

You will then need to add some extra channels (repositories to search for packages), in order to be able to create the environments:
```bash
conda config --add channels conda-forge bioconda
```

It is also recommended to install Mamba, a conda 'drop-in' replacement, before creating the environments to speed things up:
```bash
conda install -y mamba
```


At this point you might normally start creating your environments with the software packages and dependencies that you intend to use. However, Metawrap and METABOLIC are a bit complicated and rely on setting some paths to databases as well as running some setup scripts, which would take too long to do here. A separate Markdown document (@ name.md) is provided for those who wish to do this themselves at a later date.

**We will be using some pre-configured environments in this workshop that have everything you need to run these tools.**

## Check the environments are set up correctly

Although we are providing pre-configured environments, we should check to make sure they are working properly for you and can find the databases required.

### Check Metawrap environment

A Metawrap environment is provided for this training at `/uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/Metawrap-v1.3.2`. You can activate it using:
```bash
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/Metawrap-v1.3.2
```

You should see the environment (in brackets) before your username in the command prompt change as below:
![(base) changes to (Metawrap-v1.3.2) before the command prompt](images/Metawrap-env.png)

Now we need to check that the paths to the databases are configured properly. Verify that the output of:
```bash
echo $BLASTDB
```
matches:
> /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/databases/NCBI_nt

and that the output of:
```bash
echo $TAXDUMP
```
matches:
> /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/databases/NCBI_tax

If it does, then Metawrap is ready to be used.

**Deactivate the environment using:**
```bash
conda deactivate
```


### Check METABOLIC environment

A Metawrap environment is provided for this training at `/uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/METABOLIC-v4.0`. You can activate it using:
```bash
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/METABOLIC-v4.0
```

You should see the environment (in brackets) before your username in the command prompt change as below:
![(base) changes to (METABOLIC-v4.0) before the command prompt](images/METABOLIC-env.png)

Now we need to check that the paths to the databases are configured properly. Verify that the output of:
```bash
echo $GTDBTK_DATA_PATH
```
matches:
> /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/databases/GTDB-Tk_r220

and the output of:
```bash
echo $GTDB_DATA_PATH
```
matches:
> /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/databases/GTDB-Tk_r220

If it does, then METABOLIC is ready to be used.

**Deactivate the environment using:**
```bash
conda deactivate
```

## Clone this repository to your own space on Maxwell

This will give you environment files and the scripts needed to run the analysis (which will be uploaded AFTER the workshop, to force you to write them yourself for now.)

Navigate to somewhere sensible
```bash
cd /uoa/scratch/users/your-username
```

```bash
git clone https://github.com/biomoff/UoA_Metagenomics_Workshop
```

@ WHAT DOES THIS NOW PROVIDE?

---

# Get the Raw Data

A small metagenomics dataset has been prepared for this workshop. Navigate to the workshop folder
```bash
cd /uoa/scratch/users/your-username/UoA_Metagenomics_Workshop
```

and copy over the raw data:
```bash
cp -r /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/0_RAW-READS 0_RAW-READS
```

You should now have a directory called `0_RAW-READS` containing some metagenomics sequence files:

> A_1.fastq  
> A_2.fastq  
> B_1.fastq  
> B_2.fastq  


Verify that what you have matches the above:
```bash
ls -l 0_RAW-READS
```

*If it matches, you can now begin the pipeline, starting with some quality control*

---

# Quality Control using Metawrap's Read QC module

Before we can proceed with any sort of analysis, we need to do some quality control on the data we have.

We will be using the `read_qc` module of Metawrap to trim the reads. Often you might want to remove reads from a host (e.g. human), but we will be skipping that step here for simplicity as it is not required in our dataset. If you need to do it on your own data, you can find the information [here](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md).


## Making the QC.sh submission script
We need to create the `QC.sh` script to run read trimming on our two samples. Start by creating and opening this file in the `Scripts` subdirectory:
```bash
module load nano
nano Scripts/QC.sh
```

You should now be inside an open empty file, which we will gradually add to in order to build up our submission script for read QC.

### Adding Slurm scheduler parameters

We need to tell the Slurm scheduler how much resource we want to allocate to the job. The following parameters have been optimised for this dataset. Add them to your open `QC.sh` file by copying and pasting:
```bash
#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=



```
*note: add in your email address after `mail-user=` in order to get email updates for the job*

### Adding some description of the version the script was written for

It is also good practice to note which version of the package you are using the script was written for, in the event that a newer version becomes available. You can add something like the below which will be useful if you use the script later or send it to somebody else to use. Copy and paste it into your open file:
```bash
### This was written for metawrap version 1.3.2 ###


```

### Adding code for loading miniconda and activating the correct
Next, we need to make sure that our job is using the environment that we made, so that it has access to Metawrap. Add the below to your open file:
```bash
## Source own miniconda3 installation

source /uoa/home/your-username/miniconda3/etc/profile.d/conda.sh
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/Metawrap-v1.3.2


```
*note: change 'your-username' to your actual user name. If you fail to do this, conda will not work.*

### Adding code to set the input and output directories
It is also useful to be able to see at the top of the file what the input and output directories that the script will look for files in / output files to. This means that if you want to change the name of the output directory you only need to do it once at the top by changing the variable, and it will automatically apply to the whole script. Add the below to your open file:
```bash
## Set directory and file name variables

InputDir=0_RAW-READS
OutputDir=1_READ-QC


```
*note: avoid changing the names here, as the rest of this workshop will reference these directories specifically.*

### Add code to make the output directories

```bash
## Make output directories

mkdir -p "$OutputDir"
```


### Adding the code to run Metawrap's QC module

Finally you can add the code to use Metawrap to perform read trimming on your samples:
```bash
## QC and trim raw reads:

metawrap read_qc -1 "$InputDir"/A_1.fastq -2 "$InputDir"/A_2.fastq -t 10 -o "$OutputDir"/A --skip-bmtagger
metawrap read_qc -1 "$InputDir"/B_1.fastq -2 "$InputDir"/B_2.fastq -t 10 -o "$OutputDir"/B --skip-bmtagger


```

This will read in the forward reads (`-1`) and the reverse reads (`-2`) from your samples, use 10 cpu threads (`-t 10`) to perform the read QC, and then send the output files to your specified output directory (`-o`). We have used the `--skip-bmtagger` flag to skip removal of host sequences.

**Close and save the file by using Ctrl+x, then y, then Enter.**

## Running the QC.sh script

You can now run the script as below:
```bash
sbatch Scripts/QC.sh
```

This should take roughly 20-25 minutes to run.

## Inspecting the output of Metawrap's QC module

You should, once it finishes running, have an output directory called `1_READ-QC` that contains two subdirectories, `A` and `B`. These will each contain 2 files `final_pure_reads_1.fastq` and `final_pure_reads_2.fastq` for the trimmed forward and reverse reads respectively, and some output directories with QC reports: `pre-QC_report` and `post-QC_report`. Verify that you have these outputs using:
```bash
ls -l 1_READ-QC/A
ls -l 1_READ-QC/B
```

Within `pre-QC_report` and `post-QC_report` you would find:
> final_pure_reads_1_fastqc.html  
> final_pure_reads_2_fastqc.html  

These QC report HTML files can be inspected by copying them over to your own machine (using WinSCP or other methods) and viewing them in a browser. We will skip that with our data for now, but below is an example of what you might expect to see:

pre-QC reads:
![pre-QC reads](images/pre-qc.png)

post-QC reads:
![post-QC reads](images/post-qc.png)

---

# Assembling the reads into Contigs using Metawrap's Assembly module


Now that we have high quality reads, we need to assemble them into contigs using Metawrap's `assembly` module.

Depending on your study, you might want to assemble the reads from all samples together, or separately. If, for example, you have multiple replicates from a single sample (e.g. a particular soil), it may be of use to coassemble all those samples together. In our case we will stick with assembling our samples individually. 

Metawrap offers 2 different assembly algorithms, [MetaSPADES](https://doi.org/10.1101/gr.213959.116) and [Megahit](https://doi.org/10.1093/bioinformatics/btv033). MetaSPADES is usually recommended for all but large datasets, but we will use Megahit here to reduce computational burden.


## Making the Assembly.sh submission script

We need to create the `Assembly.sh` script to assemble our two samples. Start by creating and opening this file in the `Scripts` subdirectory:
```bash
nano Scripts/Assembly.sh
```

You should now be inside an open empty file, which we will gradually add to in order to build up our submission script for assembly of reads into contigs.

### Adding Slurm scheduler parameters

We need to tell the Slurm scheduler how much resource we want to allocate to the job. The following parameters have been optimised for this dataset. Add them to your open `Assembly.sh` file by copying and pasting:
```bash
#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=250M
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


```
*note: add in your email address after `mail-user=` in order to get email updates for the job*

### Adding some description of the version the script was written for

It is also good practice to note which version of the package you are using the script was written for, in the event that a newer version becomes available. You can add something like the below which will be useful if you use the script later or send it to somebody else to use. Copy and paste it into your open file:
```bash
### This was written for metawrap version 1.3.2 ###


```

### Adding code for loading miniconda and activating the correct
Next, we need to make sure that our job is using the environment that we made, so that it has access to Metawrap. Add the below to your open file:
```bash
## Source own miniconda3 installation

source /uoa/home/your-username/miniconda3/etc/profile.d/conda.sh
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/Metawrap-v1.3.2


```
*note: change 'your-username' to your actual user name. If you fail to do this, conda will not work.*

### Adding code to set the input and output directories
As above, add the below to your open file:
```bash
## Set directory and file name variables

InputDir=1_READ-QC
OutputDir=2_ASSEMBLY


```
*note: avoid changing the names here, as the rest of this workshop will reference these directories specifically.*

### Add code to make the output directories

```bash
## Make output directories

mkdir -p "$OutputDir"/A
mkdir -p "$OutputDir"/B


```

### Add code to assemble the reads into contigs using Metawrap's Assembly module

We want to run the assembly module individually on each sample, so we will do that by adding the following:
```bash
## Assemble contigs for sample A using metawrap assembly:

metawrap assembly -1 "$InputDir"/A/final_pure_reads_1.fastq -2 "$InputDir"/A/final_pure_reads_2.fastq -m 8 -t 16 --megahit -o "$OutputDir"/A

## Assemble contigs for sample B using metawrap assembly:

metawrap assembly -1 "$InputDir"/B/final_pure_reads_1.fastq -2 "$InputDir"/B/final_pure_reads_2.fastq -m 8 -t 16 --megahit -o "$OutputDir"/B


```

This will take in the trimmed forward reads `-1` and the trimmed reverse reads `-2`, tell the assembler how much memory we have allocated `-m 4` (4 GB) and how many CPU threads to use `-t 16`, to use the megahit assembly algorithm `--megahit`, and where to send the output `-o`. There is the option to change the minimum length of assembled contigs to output, `-l`, but we will leave it at the default of 1000 bp.

**Close and save the file by using Ctrl+x, then y, then Enter.**


## Running the Assembly.sh script

You can now run the script as below:
```bash
sbatch Scripts/Assembly.sh
```

It should take about 25-30 minutes to run.


## Inspecting the output of Metawrap's assembly module

Once `Assembly.sh` has finished running, you should have a subdirectory called `2_ASSEMBLY` that contains subdirectories `A` and `B`. Each of these will contain 2 output files and 2 directories containing intermediate files:

> QUAST_out  
> assembly_report.html  
> final_assembly.fasta  
> megahit  

`megahit` and `QUAST_out` contain intermediate files. `final_assembly.fasta` is the final assembly file containing the contigs and `assembly_report.html` is the assembly report generated by QUAST.

Verify that you have these files for both samples A and B by running:
```bash
ls -l 2_ASSEMBLY/A
ls -l 2_ASSEMBLY/B
```

The full assembly report can be inspected by transferring it to your own machine (using WinSCP or other methods) and viewing it in a browser. The report will contain some statistics about the assembly as well as some handy plots of cumulative assembly size by number of contigs, and the size of the contigs in descending order, as well as GC content:
![Assembly Report from QUAST](images/assembly-report.png)

However, we will just quickly inspect the output here on the command line. Copy and paste the below into the command line:
```bash
tail -n 9 2_ASSEMBLY/A/QUAST_out/report.txt
tail -n 9 2_ASSEMBLY/B/QUAST_out/report.txt
```

You should get an output that looks like this for sample A:
> \# contigs                   345  
> Largest contig              608195  
> Total length                17514579  
> GC (%)                      65.77  
> N50                         100354  
> N75                         56802  
> L50                         47  
> L75                         105  
> \# N's per 100 kbp           0.00  

and like this for sample B:
> \# contigs                   355  
> Largest contig              393717  
> Total length                16614804  
> GC (%)                      63.61  
> N50                         93956  
> N75                         54588  
> L50                         55  
> L75                         113  
> \# N's per 100 kbp           0.00  

**Now that we have assembled our reads into contigs, it is time to bin the contigs into different bins to make Metagenome Assembled Genomes (MAGs)**

---

# Binning into MAGs using Metawrap's Binning module

Binning is a crucial step in assembling MAGs from metagenomics data, and involves the placing of contigs into different bins based on taxonomy-informed or taxonomic indepedent methods. Taxonomic indepedent methods use GC content, K-mer frequencies, read depth, and co-variation of abundance across different samples to assign contigs to bins. Here we will be using [Maxbin2](https://doi.org/10.1093/bioinformatics/btv638), which uses tetranucleotide frequences, and [MetaBAT2](https://doi.org/10.7717/peerj.7359), which uses a graph-based approach on contig similarity.

*Note: It is useful to remember that these bins often represent a collection of closely related taxa. 1 MAG != 1 genome from 1 organism.*


Metawrap uses multiple different binning algorithms and combines their output to yield more accurate bins:
![Metawrap has multiple binning options that it consolidates and refines](images/binning.png)


**We will now write the submission script to run Metawrap's Binning module.**


## Making the Binning.sh submission script

We need to create the `Binning.sh` script to assemble our two samples. Start by creating and opening this file in the `Scripts` subdirectory:
```bash
nano Scripts/Binning.sh
```

You should now be inside an open empty file, which we will gradually add to in order to build up our submission script for assigning contigs to different bins.

### Adding Slurm scheduler parameters

We need to tell the Slurm scheduler how much resource we want to allocate to the job. The following parameters have been optimised for this dataset. Add them to your open `.sh` file by copying and pasting:
```bash
#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


```
*note: add in your email address after `mail-user=` in order to get email updates for the job*

### Adding some description of the version the script was written for

It is also good practice to note which version of the package you are using the script was written for, in the event that a newer version becomes available. You can add something like the below which will be useful if you use the script later or send it to somebody else to use. Copy and paste it into your open file:
```bash
### This was written for metawrap version 1.3.2 ###

```

### Adding code for loading miniconda and activating the correct
Next, we need to make sure that our job is using the environment that we made, so that it has access to Metawrap. Add the below to your open file:
```bash
## Source own miniconda3 installation

source /uoa/home/your-username/miniconda3/etc/profile.d/conda.sh
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/Metawrap-v1.3.2


```
*note: change 'your-username' to your actual user name. If you fail to do this, conda will not work.*

### Adding code to set the input and output directories
As above, add the below to your open file:
```bash
## Set directory and file name variables

Reads=1_READ-QC
Assembly=2_ASSEMBLY
OutputDir=3_INITIAL-BINNING


```
*note: avoid changing the names here, as the rest of this workshop will reference these directories specifically.*

### Add code to make the output directories

```bash
## Make output directories

mkdir -p "$OutputDir"/A
mkdir -p "$OutputDir"/B


```

### Add code to assign the contigs to bins using Metawrap's binning module

We want to run the assembly module individually on each sample, so we will do that by adding the following:
```bash
## Bin contigs using Metawrap binning for sample A

metawrap binning -o "$OutputDir"/A -t 4 -a "$Assembly"/A/final_assembly.fasta --metabat2 --maxbin2 --universal "$Reads"/A/*_1.fastq "$Reads"/A/*_2.fastq

## Bin contigs using Metawrap binning for sample B

metawrap binning -o "$OutputDir"/B -t 4 -a "$Assembly"/B/final_assembly.fasta --metabat2 --maxbin2 --universal "$Reads"/B/*_1.fastq "$Reads"/B/*_2.fastq


```

Here we are using Metawrap's `binning` module to bin contigs, specifying the number of threads (`-t 4`), the location of the assembly file (`-a`), which binning algorithms to use (`--metabat2 --maxbin2`), and to use universal marker genes rather than just bacterial marker genes (`--universal`) which improves archaeal binning. We then specify that path to the raw reads as positional arguments at the end. There are other options that can be set, but we have left those with the defaults. For your own data, you can check them using `metawrap binning --help`.

**Close and save the file by using Ctrl+x, then y, then Enter.**

## Running the Binning.sh script

You can now run the script as below:
```bash
sbatch Scripts/Binning.sh
```

It should take about 15-20 minutes to run.


## Inspecting the output of Metawrap's Binning module

You should now have a directory called `3_INITIAL-BINNING` that contains subdirectories `A` and `B`, which will each contain a subdirectory for each of the binning methods and a subdirectory with intermediate files:
> maxbin2_bins  
> metabat2_bins  
> work_files

Verify you have these using:
```bash
ls -l 3_INITIAL-BINNING/A
ls -l 3_INITIAL-BINNING/B
```

Within each of these you should have the fasta format files of the various bins contigs have been assigned to, e.g. for sample A:
```bash
ls -l 3_INITIAL-BINNING/A/maxbin2_bins
```
> bin.0.fa  
> bin.1.fa  
> bin.2.fa

```bash
ls -l 3_INITIAL-BINNING/A/metabat2_bins
```
> bin.1.fa  
> bin.2.fa  
> bin.3.fa  
> bin.unbinned.fa

**We will refine these in the next step to consolidate the results of the 2 different binners.**

---

# Refining bins using Metawrap's Bin Refinement module

**As we mentioned above, Metawrap uses multiple binners and then combines their output to create consolidated bins that are better than those produced by any single algorithm. We will use Metawrap's `bin_refinement` module to do this.**

## Making the .sh submission script

We need to create the `Bin_refinement.sh` script to assemble our two samples. Start by creating and opening this file in the `Scripts` subdirectory:
```bash
nano Scripts/Bin_refinement.sh
```

You should now be inside an open empty file, which we will gradually add to in order to build up our submission script for the refinement of the various bin sets.

### Adding Slurm scheduler parameters

We need to tell the Slurm scheduler how much resource we want to allocate to the job. The following parameters have been optimised for this dataset. Add them to your open `Bin_refinement.sh` file by copying and pasting:
```bash
#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


```
*note: add in your email address after `mail-user=` in order to get email updates for the job*

### Adding some description of the version the script was written for

It is also good practice to note which version of the package you are using the script was written for, in the event that a newer version becomes available. You can add something like the below which will be useful if you use the script later or send it to somebody else to use. Copy and paste it into your open file:
```bash
### This was written for metawrap version 1.3.2 ###

```

### Adding code for loading miniconda and activating the correct
Next, we need to make sure that our job is using the environment that we made, so that it has access to Metawrap. Add the below to your open file:
```bash
## Source own miniconda3 installation

source /uoa/home/your-username/miniconda3/etc/profile.d/conda.sh
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/Metawrap-v1.3.2


```
*note: change 'your-username' to your actual user name. If you fail to do this, conda will not work.*

### Adding code to set the input and output directories

As above, add the below to your open file:
```bash
## Set directory and file name variables

Reads=1_READ-QC
Bins=3_INITIAL-BINNING
OutputDir=4_BIN-REFINEMENT


```
*note: avoid changing the names here, as the rest of this workshop will reference these directories specifically.*

### Add code to make the output directories

```bash
## Make output directories

mkdir -p "$OutputDir"/A
mkdir -p "$OutputDir"/B


```

### Add code to refine the bins from different methods using Metawrap's bin_refinement module

We want to run the assembly module individually on each sample, so we will do that by adding the following:
```bash
## Run Metawrap's bin refinement module for sample A

metawrap bin_refinement -o "$OutputDir"/A -t 4 -A "$Bins"/A/metabat2_bins/ -B "$Bins"/A/maxbin2_bins/ -c 50 -x 10

## Run Metawrap's bin refinement module for sample B

metawrap bin_refinement -o "$OutputDir"/B -t 4 -A "$Bins"/B/metabat2_bins/ -B "$Bins"/B/maxbin2_bins/ -c 50 -x 10


```

Here, we are using Metawrap's `bin_refinement` module to take in the multiple bin sets and consolidate them into final improved bins. Here we are only using MetaBAT2 (`-A`) and Maxbin2 (`-B`) bin sets to consolidate. If we had run binning with CONCOCT, we could also supply that using `-C <path-to-concoct-bins>`. It is actually possible to supply bins produced using anything binning algorithm you like. Consult the [documentation](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md) if this is something you want to do with your own data. We have also specified the number of threads to use (`-t 4`), and where to send the output (`-o`). Finally, we have used the `-c` and `-x` flags to set thresholds for completeness and contamination, which can be tailored for your own data depending on your needs. This uses the proportion of expected single copy genes that are present (completeness), and the proportion of expected single copy genes that are present more than once (contamination) in each bin. 


**Close and save the file by using Ctrl+x, then y, then Enter.**

## Running the Bin_refinement.sh script

You can now run the script as below:
```bash
sbatch Scripts/Bin_refinement.sh
```

It should take about 1 hour and 5 minutes to run.


## Inspecting the output of Metawrap's bin refinement module

This step will produce various outputs and intermediate files in the `4_Bin-refinement` directory:

> figures  
> maxbin2_bins  
> maxbin2_bins.contigs  
> maxbin2_bins.stats  
> metabat2_bins  
> metabat2_bins.contigs  
> metabat2_bins.stats  
> metawrap_50_10_bins  
> metawrap_50_10_bins.contigs  
> metawrap_50_10_bins.stats  
> work_files

verify that what you have matches the above:
```bash
ls -l 4_BIN-REFINEMENT/A
ls -l 4_BIN-REFINEMENT/B
```


`metawrap_50_10_bins` contains the actual bins (MAGs) that were produced by the module:
> bin.1.fa  
> bin.2.fa  
> bin.3.fa

Verify that you have the same number of bins:
```bash
ls -l 4_BIN-REFINEMENT/A/metawrap_50_10_bins
ls -l 4_BIN-REFINEMENT/B/metawrap_50_10_bins
```

You can also check `metawrap_50_10_bins.stats` to see the statistics and taxonomy of the bins:

```bash
cat 4_BIN-REFINEMENT/A/metawrap_50_10_bins.stats
```
> bin     completeness    contamination   GC      lineage N50     size    binner  
> bin.1   99.57   0.0     0.722   Streptomycetaceae       80779   8364312 binsA  
> bin.3   99.45   0.0     0.554   Cyanobacteria   369767  2665800 binsA  
> bin.2   98.82   1.014   0.616   Pseudomonas     117668  5996447 binsAB  

```bash
cat 4_BIN-REFINEMENT/B/metawrap_50_10_bins.stats
```
> bin     completeness    contamination   GC      lineage N50     size    binner  
> bin.2   100.0   0.0     0.310   Euryarchaeota   115106  1846511 binsB  
> bin.1   99.57   0.0     0.722   Streptomycetaceae       80779   8365551 binsA  
> bin.3   98.82   1.044   0.616   Pseudomonas     119255  6026060 binsA  

*Note: this is a synthetic metagenome constructed from multiple individually sequenced genomes, so the completeness and contamination will be far better than you are likely to generate working on your own data. Do not be concerned if your own results do not look this good.*


**Now that we consolidated the multiple bin sets from different binners, and gotten an idea of the different taxa / genomes in our samples, we can go ahead and reassemble the reads on a per-bin basis to improve them further.**

---

# Reassembling bins using Metawrap's Reassemble_bins module

The consolidated bin set can be further improved in a lot of cases by reassembling the bins by collecting the reads belonging to each bin, then reassembling them separately. We will use Metawrap's `reassemble_bins` module to accomplish this.

*Note: in the event that reassembling the bins does not improve the result, metawrap defaults to the original bin.*


## Making the Bin_reassembly.sh submission script

We need to create the `Bin_reassembly.sh` script to assemble our two samples. Start by creating and opening this file in the `Scripts` subdirectory:
```bash
nano Scripts/Bin_reassembly.sh
```

You should now be inside an open empty file, which we will gradually add to in order to build up our submission script for reassembling bins.

### Adding Slurm scheduler parameters

We need to tell the Slurm scheduler how much resource we want to allocate to the job. The following parameters have been optimised for this dataset. Add them to your open `Bin_reassembly.sh` file by copying and pasting:
```bash
#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


```
*note: add in your email address after `mail-user=` in order to get email updates for the job*

### Adding some description of the version the script was written for

It is also good practice to note which version of the package you are using the script was written for, in the event that a newer version becomes available. You can add something like the below which will be useful if you use the script later or send it to somebody else to use. Copy and paste it into your open file:
```bash
### This was written for metawrap version 1.3.2 ###

```

### Adding code for loading miniconda and activating the correct
Next, we need to make sure that our job is using the environment that we made, so that it has access to Metawrap. Add the below to your open file:
```bash
## Source own miniconda3 installation

source /uoa/home/your-username/miniconda3/etc/profile.d/conda.sh
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/Metawrap-v1.3.2


```
*note: change 'your-username' to your actual user name. If you fail to do this, conda will not work.*

### Adding code to set the input and output directories
As above, add the below to your open file:
```bash
## Set directory and file name variables

Reads=1_READ-QC
Bins=4_BIN-REFINEMENT
OutputDir=5_BIN-REASSEMBLY


```
*note: avoid changing the names here, as the rest of this workshop will reference these directories specifically.*

### Add code to make the output directories

```bash
## Make output directories

mkdir -p "$OutputDir"/A
mkdir -p "$OutputDir"/B


```

### Add code to reassemble the reads into contigs using Metawrap's reassemble_bins module

We want to run the assembly module individually on each sample, so we will do that by adding the following:
```bash
## Run Metawrap's bin reassembly module for sample A

metawrap reassemble_bins -o "$OutputDir"/A -1 "$Reads"/A/final_pure_reads_1.fastq -2 "$Reads"/A/final_pure_reads_2.fastq -t 16 -m 64 -c 50 -x 10 -b "$Bins"/A/metawrap_50_10_bins

## Run Metawrap's bin reassembly module for sample B

metawrap reassemble_bins -o "$OutputDir"/B -1 "$Reads"/B/final_pure_reads_1.fastq -2 "$Reads"/B/final_pure_reads_2.fastq -t 16 -m 64 -c 50 -x 10 -b "$Bins"/B/metawrap_50_10_bins


```

Here we are using Metawrap's `reassemble_bins` module to reassemble the forward (`-1`) and reverse (`-2`) reads using 16 threads (`-t 16`) and 64GB memory (`-m 64`) based on the refined bins (`-b`). We are also specifying the bin completion (`-c 50`) and contamination (`-x 10`) thresholds the same as for the original binning, and specifying where the output should go (`-o`).

**Close and save the file by using Ctrl+x, then y, then Enter.**

## Running the Bin_reassembly.sh script

You can now run the script as below:
```bash
sbatch Scripts/Bin_reassembly.sh
```

It will likely take too long to run within the session, but the results are provided at `<location>`. Copy them to your own working directory with:
```bash
cp <some-location>/BIN_REASSEMBLY
```

## Inspecting the output of Metawrap's bin_reassembly module

You should now have directory called `BIN_REASSEMBLY` that contains the subdirectories `A` and `B`. Within each of these you will have the following:
> original_bins  
> original_bins.stats  
> reassembled_bins  
> reassembled_bins.checkm  
> reassembled_bins.png  
> reassembled_bins.stats  
> reassembly_results.eps  
> reassembly_results.png  
> total.faa  
> work_files  

Within `reassembled_bins.stats` you will find which version of each bin was used and the associated statistics:
```bash
cat BIN_REASSEMBLY/A/reassembled_bins.stats  
```
> bin     completeness    contamination   GC      lineage N50     size  
> bin.1.strict    99.72   0.679   0.553   Cyanobacteria   424257  2819757  
> bin.2.orig      99.57   0.0     0.721   Streptomycetaceae       84144   8530757  
> bin.3.orig      98.82   1.014   0.616   Pseudomonas     152991  5986269

```bash
cat BIN_REASSEMBLY/B/reassembled_bins.stats  
```
> bin     completeness    contamination   GC      lineage N50     size  
> bin.1.orig      98.82   1.884   0.616   Pseudomonas     127079  6043742  
> bin.2.orig      99.57   0.0     0.721   Streptomycetaceae       84144   8531397  
> bin.3.permissive        100.0   0.0     0.311   Euryarchaeota   194020  1856981

You can see that for bin 1 in sample A, and bin 3 in sample B, reassembly improved the quality of the bins. You can find the final bin files in fasta format (and some other formats), in the `reassembled_bins` subdirectory:
```bash
ls -l BIN_REASSEMBLY/A/reassembled_bins/*.fasta
```
> BIN_REASSEMBLY/A/reassembled_bins/bin.1.strict.fasta  
> BIN_REASSEMBLY/A/reassembled_bins/bin.2.orig.fasta  
> BIN_REASSEMBLY/A/reassembled_bins/bin.3.orig.fasta

```bash
ls -l BIN_REASSEMBLY/B/reassembled_bins/*.fasta
```
> BIN_REASSEMBLY/B/reassembled_bins/bin.1.orig.fasta  
> BIN_REASSEMBLY/B/reassembled_bins/bin.2.orig.fasta  
> BIN_REASSEMBLY/B/reassembled_bins/bin.3.permissive.fasta  

**These are our MAGs that we can now use for downstream analysis.**

--- 

# Functionally annotating MAGs using METABOLIC

[METABOLIC](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8) allows reconstruction of cell metabolism using genomes (including MAGs) as input. This can be done on both a community and individual genome level. Here we will be using METABOLIC in 'community' mode (`METABOLIC-C`), which should allow us to visualise the differences in function between our communities.

However, first we need to change some file suffixes from `.fa` to `.fasta` so that METABOLIC can find our sequence files. We could just do this manually for each file with `mv file.fa file.fasta`, but a script has been provided to do it in a more automated way. Run:
```bash
bash Scripts/fa_to_fasta.sh
```

You will now have the right files for METABOLIC, in a new directory called `MAGs` and we can move onto the next step.


## Making the METABOLIC.sh submission script

We need to create the `METABOLIC.sh` script to assemble our two samples. Start by creating and opening this file in the `Scripts` subdirectory:
```bash
nano Scripts/METABOLIC.sh
```

You should now be inside an open empty file, which we will gradually add to in order to build up our submission script for reassembling bins.

### Adding Slurm scheduler parameters

We need to tell the Slurm scheduler how much resource we want to allocate to the job. The following parameters have been optimised for this dataset. Add them to your open `METABOLIC.sh` file by copying and pasting:
```bash
#!/bin/bash

#SBATCH --partition=uoa-compute
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
```
*note: add in your email address after `mail-user=` in order to get email updates for the job*

### Adding some description of the version the script was written for

It is also good practice to note which version of the package you are using the script was written for, in the event that a newer version becomes available. You can add something like the below which will be useful if you use the script later or send it to somebody else to use. Copy and paste it into your open file:
```bash
### This was written for METABOLIC v4.0 ###
```

### Adding code for loading miniconda and activating the correct
Next, we need to make sure that our job is using the environment that we made, so that it has access to Metawrap. Add the below to your open file:
```bash
## Source own miniconda3 installation

source /uoa/home/your-username/miniconda3/etc/profile.d/conda.sh
conda activate /uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/envs/METABOLIC-v4.0
```
*note: change 'your-username' to your actual user name. If you fail to do this, conda will not work.*

### Adding code to set the input and output directories
As above, add the below to your open file:
```bash
## Set directory and file name variables

Metabolic=/uoa/scratch/shared/Soil_Microbiology_Group/Training/Metagenomics/METABOLIC ### WILL NEED TO CHANGE 
GenomesA=MAGs/A
GenomesB=MAGs/B
OutputDir=METABOLIC
```
*note: avoid changing the names here, as the rest of this workshop will reference these directories specifically.*

### Add code to make the output directories

```bash
## Make output directories

mkdir -p "$OutputDir"/A
mkdir -p "$OutputDir"/B
```

### Add code to specify paths to metagenomics reads

METABOLIC requires you to specify a file that gives the path to the input reads `-r`, rather than being able to pass those paths into the METABOLIC call directly. To simplify things we have written some short commands that automatically generate these files (which are then deleted once METABOLIC has run). Copy and paste the below into your open `METABOLIC.sh` file:
```bash
## Make file specifying path to metagenomic reads

echo "#Read pairs:" > A.txt
echo "READ_QC/A/final_pure_reads_1.fastq,READ_QC/A/final_pure_reads_2.fastq" >> A.txt
echo "#Read pairs:" > B.txt
echo "READ_QC/B/final_pure_reads_1.fastq,READ_QC/B/final_pure_reads_2.fastq" >> B.txt
```

### Add code to functionally annotate your genomes and samples using METABOLIC-C

```bash 
## Run METABOLIC-C on genomes for sample A
perl "$Metabolic"/METABOLIC-C.pl -t 40 -m-cutoff 0.75 -in-gn "$GenomesA" -r A.txt -kofam-db full -o "$OutputDir"/A

## Run METABOLIC-C on genomes for sample B
perl "$Metabolic"/METABOLIC-C.pl -t 40 -m-cutoff 0.75 -in-gn "$GenomesB" -r B.txt -kofam-db full -o "$OutputDir"/B
```

Here we are using METABOLIC in community mode `METABOLIC-C`, with 40 threads (`-t 40`). ADD MORE DETAIL


### Finally, clean up those files we had to make specifying the path to the reads
```bash
## Clean up
rm A.txt B.txt
```

**Close and save the file by using Ctrl+x, then y, then Enter.**

## Running the METABOLIC.sh script

You can now run the script as below:
```bash
sbatch Scripts/METABOLIC.sh
```

It should take roughly 1 hour to run.


## Inspecting the output of METABOLIC

METABOLIC provides quite a lot of output, but we will look at a subset of that here for simplicity. If you want to learn more about all the outputs it can provide, have a look at the [original publication](https://doi.org/10.1186/s40168-021-01213-8) or the [Github page](https://github.com/AnantharamanLab/METABOLIC/tree/master?tab=readme-ov-file).

You should have the following files and directories in the `METABOLIC/A` and `METABOLIC/B` output directories:
> All_gene_collections_mapped.depth.txt  
> Each_HMM_Amino_Acid_Sequence  
> KEGG_identifier_result  
> METABOLIC_Figures  
> METABOLIC_Figures_Input  
> METABOLIC_log.log  
> METABOLIC_result.xlsx  
> METABOLIC_result_each_spreadsheet  
> METABOLIC_run.log  
> MW-score_result  
> intermediate_files


 We will be looking at the following outputs:
> METABOLIC_Figures/Nutrient_Cycling_Diagrams/draw_biogeochem_cycles/draw_carbon_cycle_total.pdf  
>