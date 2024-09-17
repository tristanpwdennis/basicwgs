
# basic wgs

This repo contains the barebones pipeline for read data cleaning, mapping and qc.
Ultimately you will end up with some trimmed read files, some bam alignments, and a bunch of qc data from fastqc, qualimap, flagstat wrapped into a multiqc report

There are detailed instructions below about wrapping this into a nextflow workflow, if you wish to use nextflow.

The 'cheaper' version:

* Import the conda environment:
`conda conda env create -f <envfilename>.yml`

* Download reference genome of your choice into the `ref` directory
* Index it with `bash index.sh` in the `bin` directory
* Run `bash pipeline.sh`

### Getting started with nextflow

* First, follow the instructions below to create the docker container and run in nextflow for HPC or on local for small data.

This project uses Docker to manage all the dependencies, and nextflow to run the analysis. To get started, make sure you have docker installed. Installation instructions by platform are here:
https://docs.docker.com/engine/install/ . Once you're finished, fire up the terminal and doublecheck with ```docker -v```
Then, install nextflow: https://www.nextflow.io/docs/latest/getstarted.html . And check in the terminal again to make sure it runs ok.

### Setup and Docker installations

Clone the repo
```
git clone https://github.com/tristanpwdennis/bactocap.git
```
Enter the repo
```cd bactocap```

Now we need to build the custom Docker image for this project, and also download the GATK Docker image.
This command will build the dennistpw/align Docker image. This will take a few minutes.
```
DOCKER_BUILDKIT=1 docker build -t dennistpw/align --no-cache . 
```
Next we need to pull the GATK docker image
```
docker pull broadinstitute/gatk
```
Now let's check to make sure both of the images are ok
```
docker images
```
You should see the gatk and align repos are in the list.

### Data
When the project is further advanced, the read data for the bactocap project will be available on sra, and I will include fastq-dump commands in the workflow that will enable download of the data. Right now, however, we will have to make do the cheapo way, so please put some trimmed read files of your own into the raw_reads subdirectory, located ```datasets/<dataset>/raw_reads``` - choose whichever is appropriate for your project.

### Running the workflow
It's as simple as running 
```
nextflow run main.nf --help
```
Nextflow caches all the steps, so you don't have to go back to square one with each reanalysis. Just add more data to the raw_reads directory, or restart if you accidentally shut off your machine with
```
nextflow run main.nf -resume
```

Note, I quite like running these scripts in screen sessions: https://linuxize.com/post/how-to-use-linux-screen/
This allows me to run the workflow, check on it periodically as it runs on the other screen, whilst I tool about doing other stuff. It also reduces the likelihood of that scenario where you accidentally close your laptop when you have a terminal session running and halt your analysis - TD


### Output

The final bam files  and most qc data will be published in the ```results``` directory in each dataset directory according to sample name
The individual fastqc and bamqc data will be published in the ``` individual_reports``` subdirectory and agglomerated in the ```multiqc_report.html``` document.
A tab delimited text file ```mapping_stats.csv``` contains the flagstat data for analysis etc

