## What is it ?

A pipeline to process scRNAseq data from GEO datasets.

## POC process of GSE135194

Dataset source is [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135194): HSCPs from HD and GATA2 deficiency patients.

## Caveats

### Dataset size

13 samples, there are potential issues during SRA dump (it could be necessary to download / dump seqfiles **again** in case of uncomplete / broken fastq, check any ERROR.log in the root directory, and start the pipeline again). Pipeline process peaks at more than 2 Tio of storage: watch out volumetric issues.

### Processing time


