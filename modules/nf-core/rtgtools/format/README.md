# RTG Tools - format

[RTG Tools GitHub link](https://github.com/RealTimeGenomics/rtg-tools)

## Overview

Converts the contents of sequence data files (FASTA/FASTQ/SAM/BAM) into the RTG Sequence Data File (SDF) format.

## Command help

`rtg format --help`

```bash
Usage: rtg format [OPTION]... -o SDF FILE+
                  [OPTION]... -o SDF -I FILE
                  [OPTION]... -o SDF -l FILE -r FILE

Converts the contents of sequence data files (FASTA/FASTQ/SAM/BAM) into the RTG Sequence Data File (SDF) format.

File Input/Output
  -f, --format=FORMAT            format of input. Allowed values are [fasta, fastq, fastq-interleaved, sam-se, sam-pe] (Default is fasta)
  -I, --input-list-file=FILE     file containing a list of input read files (1 per line)
  -l, --left=FILE                left input file for FASTA/FASTQ paired end data
  -o, --output=SDF               name of output SDF
  -p, --protein                  input is protein. If this option is not specified, then the input is assumed to consist of nucleotides
  -q, --quality-format=FORMAT    quality data encoding method used in FASTQ input files (Illumina 1.8+ uses sanger). Allowed values are [sanger, solexa,
                                 illumina] (Default is sanger)
  -r, --right=FILE               right input file for FASTA/FASTQ paired end data
      FILE+                      input sequence files. May be specified 0 or more times

Filtering
      --duster                   treat lower case residues as unknowns
      --exclude=STRING           exclude input sequences based on their name. If the input sequence contains the specified string then that sequence is
                                 excluded from the SDF. May be specified 0 or more times
      --select-read-group=STRING when formatting from SAM/BAM input, only include reads with this read group ID
      --trim-threshold=INT       trim read ends to maximise base quality above the given threshold

Utility
      --allow-duplicate-names    disable checking for duplicate sequence names
  -h, --help                     print help on command-line flag usage
      --no-names                 do not include name data in the SDF output
      --no-quality               do not include quality data in the SDF output
      --sam-rg=STRING|FILE       file containing a single valid read group SAM header line or a string in the form
                                 "@RG\tID:READGROUP1\tSM:BACT_SAMPLE\tPL:ILLUMINA"
```
