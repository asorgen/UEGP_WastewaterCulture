
# UEGP Wastewater Culturing Analysis

## Transparency and Reproducibility

This repository contains:         

 * Counts files produced by processing the sequences, see folder: `data`
 * Configuration file to run the analysis component of the project, including all Rcode: `otu-analysis.properties`
 * Results of the analysis (figures and tables), see folder: `results`
 
See `RScripts/README` for instructions on how to repeat the entire analysis on your own machine, thus reproducing all figures and tables.  A review module was added to the pipeline to facilitate readers and reviewers comparing the reproduced output the the corresponding components in the manuscript.

## Manuscript

**Title:**                
_Characterization of Environmental and Cultivable Antibiotic-Resistant Microbial Communities Associated with Wastewater Treatment_

**Authors:**               
Alicia Sorgen, James Johnson, Kevin Lambirth, Sandra M. Clinton, Molly Redmond, Anthony Fodor, and Cynthia Gibas

_Antibiotics_ 2021
www.mdpi.com/journal/antibiotics


## Abstract
Bacterial resistance to antibiotics is a growing global concern, threatening human and environ-mental health, particularly among urban populations. Wastewater treatment plants (WWTPs) are thought to be “hotspots” for antibiotic resistance dissemination. The conditions of WWTPs, in conjunction with the persistence of commonly used antibiotics, may favor the selection and trans-fer of resistance genes among bacterial populations. WWTPs provide an important ecological niche to examine the spread of antibiotic resistance. We used heterotrophic plate count methods to identify phenotypically resistant cultivable portions of these bacterial communities and charac-terized the composition of the culturable subset of these populations. Resistant taxa were more abundant in raw sewage and wastewater before the biological aeration treatment stage. While some antibiotic-resistant bacteria (ARB) were detectable downstream of treated wastewater re-lease, these organisms are not enriched relative to effluent-free upstream water, indicating effi-cient removal during treatment. Combined culture-dependent and -independent analyses re-vealed a stark difference in community composition between culturable fractions and the envi-ronmental source material, irrespective of culturing conditions. Higher proportions of the envi-ronmental populations were recovered than predicted by the widely accepted 1% culturability paradigm. These results represent baseline abundance and compositional data for ARB commu-nities for reference in future studies addressing the dissemination of antibiotic resistance associ-ated with urban wastewater treatment ecosystems.


**Funding:**
This study was funded in part by a grant from NIH/NIMHD (R01MD011504), NIH/NIEHS (P30 ES010126). ER was supported by T32ES007018. The Pregnancy, Infection, and Nutrition study was supported by Grants HD37584, HD39373, and DK61981. The General Clinic Research Center was supported by the National Institutes of Health General Clinical Research Centers program of the Division of Research Resources Grant RR00046. 

**Availability of data and material:**
All sequence data have been submitted to the National Center for Biotechnology Information (NCBI) under BioProject ID number PRJNA657079. All analyses and visualizations were generated using a Dockerized version of R Studio (Version 4.0.2) and the R scripts are available at https://github.com/asorgen/UEGP_WastewaterCulture.