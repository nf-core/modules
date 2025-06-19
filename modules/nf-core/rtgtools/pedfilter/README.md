# RTG Tools - pedfilter

[RTG Tools GitHub link](https://github.com/RealTimeGenomics/rtg-tools)

## Overview

Filter and convert a pedigree file.

## Command help

`rtg pedfilter --help`

```bash
Usage: rtg pedfilter [OPTION]... FILE

Filter and convert a pedigree file.

File Input/Output
      FILE                 the pedigree file to process, may be PED or VCF, use '-' to read from stdin

Filtering
      --keep-family=STRING keep only individuals with the specified family ID. May be specified 0 or more times, or as a comma separated list
      --keep-ids=STRING    keep only individuals with the specified ID. May be specified 0 or more times, or as a comma separated list
      --keep-primary       keep only primary individuals (those with a PED individual line / VCF sample column)
      --remove-parentage   remove all parent-child relationship information

Reporting
      --vcf                output pedigree in the form of a VCF header

Utility
  -h, --help               print help on command-line flag usage
```
