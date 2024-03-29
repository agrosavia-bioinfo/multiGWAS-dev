# MultiGWAS Table of Contents
   * [MultiGWAS Table of Contents](#multigwas-table-of-contents)
   * [Installation](#installation)
      * [General steps to install multiGWAS on a Linux system](#general-steps-to-install-multigwas-on-a-linux-system)
      * [Specific instructions to install multiGWAS on a Linux Ubuntu](#specific-instructions-to-install-multigwas-on-a-linux-ubuntu)
         * [Install external software:](#install-external-software)
         * [Install MultiGWAS tool:](#install-multigwas-tool)
   * [Executing the MultiGWAS tool](#executing-the-multigwas-tool)
      * [Using the command line interface:](#using-the-command-line-interface)
      * [Using the graphical user interface:](#using-the-graphical-user-interface)
   * [Running the examples](#running-the-examples)
   * [General usage](#general-usage)
   * [Genomic data formats](#genomic-data-formats)
      * [Genotype:](#genotype)
      * [Phenotype](#phenotype)
   * [Configuration file](#configuration-file)
   * [Considerations](#considerations)
      * [Implementation](#implementation)
      * [Number of SNPs in Manhattan and QQ plots](#number-of-snps-in-manhattan-and-qq-plots)
      * [Correction for multiple testing](#correction-for-multiple-testing)

# Installation
We describe here the [MultiGWAS](https://github.com/agrosavia-bionformatics/multiGWAS) installation process that contains the file sources and binaries to run the MultiGWAS tool. It includes:
  - MultiGWAS tool sources (R code, Java application,  and binary bash scripts)
  - MultiGWAS precompiled R libraries for Linux systems (tested on Ubuntu 18.04) 
  - Binaries and Java classes for the four GWAS packages  GWASpoly, SHEsis, TASSEL, and PLINK (r1.9 and r2.0)
  
## General steps to install multiGWAS on a Linux system
  1. Install pandoc markup converter, if not installed. 
  2. Install git tool to clone github repository, if not installed.
  3. Install Oracle Java runtime, if not installed.
  4. Install R 3.6, if not installed. 
  5. Clone the multiGWAS repository
  6. Execute the multiGWAS installer:

## Specific instructions to install multiGWAS on a Linux Ubuntu
### Install external software
The MultiGWAS tool currently runs on Linux systems (tested on Ubuntu Linux 18.04 LTS, x86_64 GNU / Linux), and requires the following software to be installed:
  - Install git tool and pandoc markup converter: 
```
    sudo apt install git pandoc
```
  - R 3.6 or higher. If not installed see https://cran.r-project.org/bin/linux/ubuntu/README.html or Open a Linux console and enter the following instructions for Ubuntu 18.04 (bionic):
```
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
	sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
	sudo apt-get update
	sudo apt-get install r-base-dev
```
  - Oracle Java JRE 8.x or higher. If not installed see https://ubuntuhandbook.org/index.php/2020/09/install-oracle-java-15-ubuntu-20-04-18-04-16-04/ or Open a Linux console (or terminal) and enter the following instructions for Ubuntu 18.04 (bionic):
```
	sudo add-apt-repository ppa:linuxuprising/java
	sudo apt update
	sudo apt install oracle-java15-installer
	java -version
``` 

### Install MultiGWAS tool
Open a Linux console (or terminal) and enter the following instructions:
```
# 1. Clone the MultiGWAS repository:
    git clone https://github.com/agrosavia-bioinformatics/multiGWAS.git
# 2. Change to the multiGWAS dir:
    cd multiGWAS
# 3. Execute the installer:
    . INSTALL.SH         # dot space INSTALL.SH
# 4. The MultiGWAS tool is ready to use.
```

# Executing the MultiGWAS tool
MultiGWAS can be executed by either a command line interface (CLI) developed in R or a graphical user interface (GUI) developed in Java.
## Using the command line interface
  - Open a Linux console
  - Create a new folder (e.g. test)
  - Create a new configuration file (e.g. test.config, see below)  and save it into the created folder.
  - Run the CLI Multigwas program followed by the name of the configuration file:
      ```
      multiGWAS test.config
      ```
  - Open a file browser and view the results saved into the created folder (e.g. test/out-test/report/). They include:
     - Full html report (file multiGWAS-report.html)
     - Original graphics and tables generated by MultiGWAS.
     - Preprocessed tables with the GWAS results from the four GWAS packages
     
## Using the graphical user interface:
   - Open a linux console and simply enter
       ```
       jmultiGWAS
       ```
    - The GUI application is easy and straighforward. It includes four views:
       - Inputs view:  to create the configuration file and start the execution of MultiGWAS
       - Outputs view: to view the logs/messages from the current execution
       - Results view: to view a preview of the HTML report.
       - Files view:   to explore and open the results, including the original tables, graphics, and the full HTML report.

# Running the examples
The multiGWAS folder contains an "example" subfolder with genomic data ("example-genotype.tbl" and "example-phenotype.tbl" files) and two ready-to-use configuration files ("full.config" and "naive.config" files for **Naive** and **Full** GWAS analysis, respectively). Both genotype and phenotype come from the Solanaceae Coordinated Agricultural Project (SolCAP) potato diversity panel described in the paper. 

To run a Full GWAS analysis with the multiGWAS tool, follow the next steps:
 - Open a Linux terminal
 - Change to the "multiGWAS" folder
 - Change to the "example" subfolder
 - Execute the **multiGWAS tool** using as argument the configuration file "full.config":
  ```
        multiGWAS full.config 
  ```
 - An output folder will be created named as "out-XXXX" where XXXX is the prefix of the configuration filename. For the above example, an "out-full" subfolder will be created with the following files and subfolders:
    - A "multiGWAS-report.html" file in HTML format with the full report from the multiGWAS results.
    - A "report" subfolder with the resulting tables and graphics included in the previous report file.
    - An "out" subfolder that contains temporary files created by multiGWAS and the other GWAS tools.
    - A "logs" subfolder that contains the log outputs from the different tools.
    
# General usage
  - Create a new folder (e.g. "test" folder).
  - Copy the phenotype and genotype files to new folder (see data formats below)
  - Create a configuration file (e.g. test.config) and copy it to the new folder (see configuration file below)
  - Open a Linux terminal
  - Change to the new folder
  - Execute the multiGWAS tool using as argument the configuration file
```
      multiGWAS test.config
```
  - Results will be saved in the "out-test" folder

# Configuration file
The configuration file is a text file with a list of parameter names and their values separated by a colon (":"). This file is the main input for MultiGWAS and it can be created in three ways: using a general text editor, using the MultiGWAS GUI interface, or modifying an existing configuration file. 

All keywords in the parameter file may be entered in upper and/or lower case; spaces are allowed but not required around the “:” sign. 
Upper / lower case is relevant for the filenames. Blank lines and comment lines (starting with the number sign #) may be added, as well as extra text after the data on a line. The only requirement is that all the data (filenames and other strings) do not contain blanks (spaces, tab characters).

Now, we briefly describe these parametes and then we show an example of a config file.

|  Parameter name   | Description |                     
|------------------ |------------  |                   
| genotypeFile      | Genotype filename, file with the marker data (see genomic data section below) |
| genotypeFormat    | Genotype format, currently four formats: "gwaspoly", "kmatrix", "vcf", and "fitpoly" (see genomic data section below) |
| phenotypeFile     | Phenotype filename, file with the individuals and trait values (see genomic data section below) |
| mapFile           | Map file, optional file with marker information (marker, reference allele, alternate, allele, chromosome, and position (see genomic data section below) |
| significanceLevel | The genome-wide significance threshold α (commonly 0.01 or 0.05)|
| correctionMethod  | The method for multiple testing correction (”Bonferroni” or ”FDR”)|
| gwasModel         | The type of GWAS analysis (”Naive” or ”Full”)|
| nBest             | Number of top associations to be reported
| filtering         | TRUE or FALSE whether to use quality control (QC) filters or not (see below) |
| MAF               | Minor allele frequency QC filter |
| MIND              | Individual missing rate QC filter |
| GENO              | SNP missing rate QC filter |
| HWE               | Hardy-Weinberg threshold QC filter|
| tools             | Tools to be used in the analysis. Any combination of the following tools: "GWASpoly", "TASSEL", "PLINK, and SHEsis
| geneAction        | Gene-action assumed model (Marker-effect). Currently, four options: "additive", "general", "dominance", or "all" |

## Example of a configuration file
This is the contents of a typical configuration file named as "full-tetra.config":
```
genotypeFile         : example-genotype.tbl
phenotypeFile        : example-phenotype.tbl
mapFile              : 
genotypeFormat       : gwaspoly
significanceLevel    : 0.05
correctionMethod     : Bonferroni
gwasModel            : Naive
nBest                : 10
filtering            : TRUE
MAF                  : 0.01
M IND                 : 0.1
GENO                 : 0.1
HWE                  : 1e-10
tools                : GWASpoly SHEsis PLINK TASSEL
geneAction           : additive
```
## Genomic data inputs
The following parameters from configuration file are related with the type of genomic data required by MultiGWAS. Below, we show the characteristics and structure of the input files. Keep in mind that the headeer line must be present in all the file formats we show below.
### "genotypeFile"
It specifies the filename for the genotype data. (see below the accepted genotype formats).

### "genotypeFormat"
Currently, MultiGWAS accepts five genotype formats: "gwaspoly", "kmatrix", "vcf", "fitpoly", and "updog":
- ***"gwaspoly" format:*** table with comma separated values (.csv). Each row contains the marker id, the chromosome, the position in the chromosome, and the following columns correspond to the marker data for each individual codified in the "ACGT" format (e.g., AATT, CCGG, AAAT, GGCG). An example follows:
```	 
| Marker   | Chrom | Position | ACBrador | ACLPI175395 | ADGPI195204 | AdirondackBlue |
|----------|-------|----------|----------|-------------|-------------|----------------|
| c2_41437 | 0     | 805179   | AAGG     | AAAA        | AAAA        | AAGG           |
| c2_24258 | 0     | 1252430  | GGGG     | GGGG        | GGGG        | AAGG           |
| c2_21332 | 0     | 3499519  | TTCC     | CCCC        | CCCC        | TTCC           |
```

- ***"kmatrix" format:*** table with comma separated values (.csv). Each ro contains the marker id and the marker data for each individual codified in a numeric format (e.g. 0,1,2,3,4). An example follows:
 
```
| Marker   | ACBrador | ACLPI175395 | ADGPI195204 | AdirondackBlue | 
|----------|----------|-------------|-------------|----------------|
| c2_41437 | 2        | 0           | 0           | 2              |
| c2_24258 | 0        | 0           | 0           | 2              |
| c2_21332 | 2        | 4           | 4           | 2              |
```

- ***"vcf" format:*** Variang Call Format (VCF) with metadata in the first lines followed by a header line. The following lines contain genotype information of the individuals for each position. VCF marker data can be encoded as simple genotype calls (GT format field, e.g. 0/0/1/1 for tetraploids or 0/1 for diploids) or using the NGSEP custom format fields (Tello et al., 2019): ACN, ADP or BSDP. An example follows:
```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  sample01    sample02    sample03
   0   805179  c2_41437    A   G   .   .   PR  GT  0/1/1/0 0/1/1/0 0/1/0/0
   0   1252430 c2_24258    G   A   .   .   PR  GT  0/1/0/0 0/1/1/0 0/0/1/0
   0   3499519 c2_21332    T   C   .   .   PR  GT  0/1/1/0 0/1/1/1 0/1/1/0
```

- ***"fitpoly" format:*** table with tab separated values (scores file e.g. filePrefix_scores.dat) containing one line per sample for every marker that could be fitted. MultiGWAS uses only two columns from this file: "MarkerName" and "geno", with the name of the marker and the assigned genotype number, respectively. The genotype number is assigned according to the ploidy, for tetraploids from 0 to 4, and for diploids from 0 to 2. An example follows: 
```
| marker | MarkerName | SampleName     | ratio | P0    | P1    | P2    | P3    | P4    | maxgeno | maxP  | geno |
|--------|------------|----------------|-------|-------|-------|-------|-------|-------|---------|-------|------|
| 1      | c1_1       | ACBrador       | 0.932 | 0.821 | 0.407 | 0.879 | 0.537 | 0.158 | 1       | 0.363 | 1    |
| 1      | c1_1       | ACLPI175395    | 0.719 | 0.823 | 0.385 | 0.213 | 0.834 | 0.873 | 0       | 0.189 | 0    |
| 1      | c1_1       | ADGPI195204    | 0.064 | 0.09  | 0.081 | 0.477 | 0.066 | 0.845 | 0       | 0.351 | 0    |
| 1      | c1_1       | AdirondackBlue | 0.606 | 0.853 | 0.62  | 0.44  | 0.954 | 0.047 | 1       | 0.151 | 1    |
| 1      | c1_1       | AdirondackRed  | 0.785 | 0.176 | 0.079 | 0.915 | 0.634 | 0.22  | 1       | 0.259 | 1    |
```

- ***"updog" format:*** table with comma separated values containing one line per sample for every marker that could be fitted. MultiGWAS uses only two columns from this file: "snp" and "geno", with the name of the marker and the assigned genotype number, respectively. The genotype number is assigned according to the ploidy, for tetraploids from 0 to 4, and for diploids from 0 to 2. An example follows: 
```
|     snp     | ind   | ref | size | geno | postmean | maxpostprob | Pr_0 | ... | Pr_4 | logL_0  | ... | logL_4 |
|-------------|-------|-----|------|------|----------|-------------|------| ... |------|---------| ... |--------|
|PotVar0089524|P3PEM05| 113 | 143  | 3    | 2.99     | 0.99        | 0    | ... | 0    | -176.49 | ... | -39    |
|PotVar0089524|P2PEM10| 86  | 96   | 3    | 3        | 1           | 0    | ... | 0    | -172.23 | ... | -15.93 |
|PotVar0089524|P5PEM08| 122 | 142  | 3    | 3        | 1           | 0    | ... | 0    | -201.39 | ... | -25.8  |
|PotVar0089524|P7PEM09| 80  | 80   | 4    | 4        | 1           | 0    | ... | 1    | -190.47 | ... | -0.18  |
|PotVar0089524|P5PEM04| 69  | 69   | 4    | 4        | 1           | 0    | ... | 1    | -172.29 | ... | -0.16  |
```

### "phenotypeFile"
The phenotype file is formatted as a table separated by commas with the names of the variables (columns) in the first line and with two columns. The first one containing the name (or ID) of the individual, and the second the trait value. An example follows:

```
|Name           | tuber_shape |
|---------------|-------------|
|ACBrador       |    3.59     |   
|AdirondackBlue |    4.07     |
|AllBlue        |    4.73     |
|AlpineRusset   |    4.85     |
|Alturas        |    4.46     |
```

### "mapFile"
This file contains the markers information and it is an an *optional* file only needed when the genotype is numeric, as in the case of  "kmatrix" and "fitpoly" formats. The file contains a table with separated column values (.csv) with the following information for each marker: marker name, reference allele, alternate allele, chromosome, and position in the chromosome. An example follows:
```
| Markers  | Ref | Alt | Chrom | Position |
|----------|-----|-----|-------|----------|
| c2_41437 | A   | G   | 0     | 805179   |
| c2_24258 | G   | A   | 0     | 1252430  |
| c2_21332 | T   | C   | 0     | 3499519  |
| c2_21320 | T   | G   | 0     | 3810687  |
| c2_21318 | C   | T   | 0     | 3810936  |
```


# Considerations
## Implementation
Most of the code uses the R language. However, some scripts that calling the GWAS tools are writing in bash. The version of the four tools are GWASpoly 1.3 (R library), SHEsis 1.0 (binary program), PLINK 1.9 and 2.0 (binary programs), and TASSEL 5.0 (Java packages). PLINK 1.9 is used for GWAS analysis (association between SNPs and quantitative traits), and PLINK 2.0 is used to account for cryptic relatedness (estimating kinship coefficientes).

## Number of SNPs in Manhattan and QQ plots
The Manhattan and QQ plots for the different GWAS tools show a different number of markers (SNPs). Two reasons explain this pattern. First, the GWASpoly software uses four models for the marker effect (i.e., additive, general, simplex dominance, and duplex dominance). Therefore,  the plots show the SNPs four times, one for each model. Second, MultiGWAS is using scores instead of raw p-values, and scores are the -log10(p) results. So, when p-values are high, the scores have a negative value, and because the y-axes in the plot start in zero, they are not shown.

## Correction for multiple testing
MultiGWAS is using two methods for correction for multiple testing, the Bonferroni correction and adjusting the False Discovery Rate (FDR). MultiGWAS calculates the Bonferroni correction using the number of non-missing genotypes (NMISS) included in the analysis instead of the whole genotypes. Only SHEsis, PLINK, and TASSEL give the NMISS number. In contrast, GWASpoly does no show the NMISS number, but it uses it internally to calculate the corrections.

