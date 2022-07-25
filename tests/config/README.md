# Modules Test Data

Test data for modules that have been merged to [nf-core](https://github.com/nf-core/modules)
should be updated via [their SOP](https://github.com/nf-core/test-datasets/blob/modules/README.md).

Test data for Tree of Life modules live on the farm at `/lustre/scratch123/tol/resources/nextflow/test_data`.
Data should follow the usual directory structure (e.g. `genomic_data/<specimen-id>/<data-type>/`) as much
as possible, and all be registered in [test_data.config](test_data.config).

## Datasets

We include some complete data files from organisms of different genome sizes:

- small genomes: species < 100 Mb. Currently _Anthocharis cardamines_ (`ilAntCard1`), _Asterias rubens_ (`eAstRub1`), _Eimeria tenella_ (`pEimTen1`), and _Laetiporus sulphureus_ (`gfLaeSulp1`)
- medium genomes: species < 1 GB. Currently _Erannis defoliaria_ (`ilEraDefo1`), _Erithacus rubecula_ (`bEriRub1`), _Pararge aegeria_ (`ilParAegt2`)
- large genomes: species > 1 GB. Currently _Cervus elaphus_ (`mCerEla1`), _Meles meles_ (`mMelMel1`), and _Sciurus vulgaris_ (`mSciVul1`)

We also have hand-crafted test datasets, e.g. a subset of 1,000 reads, that can be more suitable for unit-testing.

## Naming guidelines

The name should include the file format and the extension.

Paired data should have identical names, with `_1` and `_2`.

## Data types

### Genomic data - `genomic_data`

The typical keys for genomic data are:

- `pacbio_bam` and `pacbio_pbi`: PacBio (HiFi) files
- `hic_cram` and `hic_crai`: Hi-C files (from any kit)
- `illumina_cram` and `illumina_crai`: Regular Illumina sequencing
- `ont_fastq_gz`: ONT FastQ file
- `tenx_i1_fastq_gz`, `tenx_r1_fastq_gz`, and `tenx_r2_fastq_gz`: 10X files from the same set (I1, R1, and R2)

Aligned data will have `aligned_` in their name, e.g. `hic_aligned_1_bam`.

### Assembly data - `assembly`

This category covers data created during the assembly pipelines.

- `hicanu_trimmed_reads_fasta_gz`: `trimmedReads.fasta` file from a HiCanu assembly directory
- `canu_contigs_fasta`: Contig-level assembly produced by Canu

### Software data

This category corresponds to other inputs and outputs for software.

- `*_longranger_mkref_targz`: archive of the output of `longranger mkref` on an assembly

## Usage

Instead of using full URLs, chain the lookup keys like this:
`params.tol_test_data['test']['dImpGla2']['genomic_data']['hic_aligned_1_bam']`

Nextflow will replace it with the actual URL it's constructed from the configuration file.

## Adding new data

The directory on disk `/lustre/scratch123/tol/resources/nextflow/test_data` can be written to by anyone in the `tolengine` group. Feel free to add your own data there, and use those in your modules. However, to get the module merged in, we require that you update the [`test_data.config`](test_data.config) file as well in your pull-request. The ToL-IT team can upload the data files to the S3 server for a final test before merging the pull-request.
