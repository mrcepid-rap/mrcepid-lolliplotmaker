# MRCEPID-Lolliplot Maker

This is a short R script that makes lolliplots with annotated variant information like the following example:

![](https://github.com/mrcepid-rap/mrcepid-lolliplotmaker/blob/main/sample_images/CHEK2_sample.png)

## Installation

### Command-line tools:

This script requires the following tools to be installed and in your `$PATH`:

* [tabix](http://www.htslib.org/doc/tabix.html) – **MUST be a version > 1.0**
* [bedtools](https://bedtools.readthedocs.io/en/latest/)

### Required R Packages:

This script requires the following R packages to be installed:

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
* [lemon](https://cran.r-project.org/web/packages/lemon/index.html)
* [bedr](https://cran.r-project.org/web/packages/bedr/index.html)
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html)

### Getting This Tool

After installing required packages/tools, do:

```
git clone https://github.com/mrcepid-rap/mrcepid-lolliplotmaker.git
cd mrcepid-lolliplotmaker/
./lolliplot_maker.R <PHENOTYPE>.bolt.markers.BOLT.stats.tsv.gz <MASK> <MAF> <ENST>
```

## Running

To run this tool, use a command line like the following:

```
./lolliplot_maker.R --markers <PHENOTYPE>.bolt.markers.BOLT.stats.tsv.gz --mask <MASK> --maf <MAF> --enst <ENST>
```

All options are described below in detail:

| short arg | long arg    | description                                                                                                                                                                                 |
|-----------|-------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -m        | --markers   | This is the per-marker output from BOLT-LMM generated by [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting).                                     |
| -a        | --mask      | The short name of the mask that you wish to plot from. See below for all possible values                                                                                                    |
| -f        | --maf       | The short name of the MAF bin that you wish to plot from. See below for all possible values                                                                                                 | 
| -e        | --enst      | An ENST Id of the gene you wish to plot. **EITHER** --enst OR --symbol must be provided                                                                                                     |
| -s        | --symbol    | The HGNC gene symbol you wish to plot. **EITHER** --enst OR --symbol must be provided                                                                                                       | 
| -p        | --phenotype | A phenotype name. This is NOT used to select data and simply adds an additional label on the x-axis of the resulting plot (see example image above).                                        |
| -o        | --output    | An output prefix to use for naming the resulting outputs. A path can be included at the beggining of this string to redirect to a different folder. This option is not required by default. |
| -h        | --help      | Print the help dialogue                                                                                                                                                                     |

### Important Notes on Inputs:

1. To properly work, the file supplied to `--markers` must still be bgzipped and tabixed indexed. If you have unzipped this file, easy commands to rezip are:

```
bgzip <PHENOTYPE>.bolt.markers.BOLT.stats.tsv
tabix -S 1 -s 2 -b 3 -e 4 <PHENOTYPE>.bolt.markers.BOLT.stats.tsv.gz
```

bgzip/tabix must be installed and in your path to work, but this should already be the case as a requirement of this tool is to have tabix.

2. --mask must be one of 9 possible variant masks you are interested in:

| Command Line Name    | Long Description                            | Query String                                                                                                  |
|----------------------|---------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| HC_PTV               | High-confidence PTVs                        | `FILTER=="PASS" & PARSED_CSQ=="PTV" & LOFTEE=="HC"`                                                           |
| PTV                  | All PTVs                                    | `FILTER=="PASS" & PARSED_CSQ=="PTV"`                                                                          |
| MISS                 | All  Missense Variants                      | `FILTER=="PASS" & PARSED_CSQ=="MISSENSE"`                                                                     |
| MISS_CADD25          | Missense Variants CADD ≥ 25                 | `FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & CADD>=25`                                                          |
| MISS_REVEL0_5        | Missense Variants REVEL ≥ 0.5               | `PARSED_CSQ=="MISSENSE" & REVEL>=0.5`                                                                         |
| MISS_REVEL0_7        | Missense Variants REVEL ≥ 0.7               | `FILTER=="PASS" & PARSED_CSQ=="MISSENSE" & REVEL>=0.7`                                                        |
| MISS_CADD25_REVEL0_7 | Missense Variants REVEN ≥ 0.7 AND CADD ≥ 25 | `PARSED_CSQ=="MISSENSE" & CADD>=25 & REVEL>=0.7`                                                              |
| SYN                  | All synonymous variants                     | `FILTER=="PASS" & PARSED_CSQ=="SYN"`                                                                          |
| DMG                  | Combined HC_PTV and MISS_CADD25             | <code>FILTER=="PASS" & ((PARSED_CSQ=="PTV" & LOFTEE=="HC") &#124; (PARSED_CSQ=="MISSENSE" & CADD>=25))</code> |

3. --maf must be one of 2 possible allele frequency bins:

| Command Line Name | Long Description               |
|-------------------|--------------------------------|
| MAF_01            | Minor Allele Frequency < 0.001 |
| AC_1              | Allele Count == 1              |

4. **EITHER** --enst or --symbol MUST be provided.
   + --enst is an ENSEMBL gene ID for the gene you are interested in. It should look like 'ENST00000000001'.
   + --symbol is an exact matching hgnc gene ID. If multiple genes with the same symbol are found, an error will be
        reported and it is recommended to re-try with an ENST as these will always be unique.

## Outputs

By default, this script generates two outputs with a file prefix like `<SYMBOL>_<MASK>_<MAF>_BOLT.*`, where `<SYMBOL>` 
is the gene SYMBOL for the gene provided with ENST/SYMBOL. This prefix can be changed with the `-o/--output` option. 

1. An image like the example given at the top of the page.
2. A tab-delimited TSV file of annotated variants that were used to build the image.

These files will be placed in the root directory of this repository.
