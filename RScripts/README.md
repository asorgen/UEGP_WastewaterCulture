# Steps to Reproduce UEGP Wastewater Culture Analysis

## 1. Ensure that BioLockJ (v1.3.13 or newer) is on local system
BioLockJ - https://biolockj-dev-team.github.io/BioLockJ/Getting-Started/

## 2. Download UEGP_WastewaterCulture directory
`git clone https://github.com/asorgen/UEGP_WastewaterCulture.git`

## 3. Set up required software

### Option A) using docker

Install docker.
Docker Desktop - https://www.docker.com/products/docker-desktop

Make sure the ` docker run hello-world ` command runs successfully.

The docker images required for this pipeline will be automatically pulled from the docker hub as needed.  The first time the pipeline runs, startup will be slow as images are downloaded. 

### Option B) not using docker

Make sure R is installed.  See https://www.r-project.org/.  These scripts were written with 4.0.2.

Make sure all required R packages are installed                                

 * data.table
 * dplyr
 * gdata
 * ggplot2
 * ggsignif
 * grDevices
 * grid
 * plotrix
 * plyr
 * scales
 * stringr
 * tidyr
 * tidyverse
 * vegan

## 4. Run BioLockJ pipeline

Move to the Analysis folder:            
`cd <path/to/UEGP_WastewaterCulture/`

To run the pipeline using **locally installed software**:                 
`biolockj out-analysis.properties`

To run the pipeline using **docker images**, add the -d argument:                                    
`biolockj -d otu-analysis.properties`


## 5. Run individual scripts

The function of each R script is provided at the beginning of each script file.
