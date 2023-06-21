# Automatically renaming files for Cellranger

See also https://github.com/nf-core/scrnaseq/issues/241 and the [Cellranger documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input).

## Motivation

Cellranger (and also Spaceranger, and probably other 10x pipelines) rely on the following input

- `--fastqs`, to point to a directory with fastq files that are named according to `{sample_name}_S{i}_L00{j}_{R1,R2}_001.fastq.gz`.
- `--sample`, with a sample name (that is the prefix of all associated fastq files)

In the easiest case, `--fastqs` points to a directory that contains all fastqs of a single sample that already follow the naming conventions. However, it's not always easy:

1.  the directory may contain fastqs for multiple samples
    - this is not a problem for cellranger, it will automatically choose the correct fastq files via the `--sample` flag as long as they follow the naming convention, **but**
    - if we stage the whole folder (or all files in that folder) using nextflow, it breaks caching in that if an additional sample (or any other file) gets added to the folder, the cache gets invalidated for _all_ samples.
2.  the sample has been sequenced across multiple flow cells
    - In this case, we need to specify multiple input folders. Cellranger allows passing a _list of directories_ e.g. `--fastqs=/path/dir1,path/dir2`
    - Staging all _files_ in these folders into a single directory using nextflow doesn't do the job, as there may be duplicate file names across the different folders.
3.  the raw sequencing data may have been downloaded from a sequence database and doesn't follow the naming convention anymore.
    - In that case we need to rename the files to follow the bcl2fastq naming convention.

## Solution

Rename files automatically to match the `{sample_name}_S{i}_L00{j}_{R1,R2}_001.fastq.gz` pattern. We assume that
files in the input channel are ordered according to `file1 R1, file1 R2, file2 R1, file2 R2, ...`.
If the files are explicitly listed in a samplesheet, we have this order anyway. If the files follow any sensible naming
convention, this order can be achieved by sorting by filename.

## Does renaming the files change the result?

No. Here's how the same dataset was tested in four different scenarios. The result is the same down
to identical md5sums of the output files.

### Prepare datasets

1. Download dataset:

```
curl https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.1.2/10k_PBMC_3p_nextgem_Chromium_X_intron/10k_PBMC_3p_nextgem_Chromium_X_intron_fastqs.tar | tar -xv
INPUT_DIR=10k_PBMC_3p_nextgem_Chromium_X_fastqs
# remove index files
rm $INPUT_DIR/*_I{1,2}_*.fastq.gz
```

2. Simulate scenario where files were generated using multiple flow cells:

```
MFC=fastq_multiflowcell
mkdir -p $MFC/fc{1,2,3,4}
for f in $INPUT_DIR/*.fastq.gz; do ; ln $f $MFC; done
for i in $(seq 4) ; do ; mv $MFC/*L00$i* $MFC/fc$i ; done
for i in $(seq 4) ; do ; rename L00$i L001 $MFC/fc$i/*; done
```

3. Simulate scenario where sample was multiplexed:

```
MPX=fastq_multiplexed
mkdir $MPX
for f in $INPUT_DIR/*.fastq.gz; do ; ln $f $MPX; done
for i in $(seq 4) ; do ; rename S1_L00$i S${i}_L001 $MPX/* ; done
```

4. Concatenate files

```
mkdir fastq_cat
cat $INPUT_DIR/*_R1_*.fastq.gz > fastq_cat/10k_PBMC_3p_nextgem_Chromium_X_S1_L001_R1_001.fastq.gz
cat $INPUT_DIR/*_R2_*.fastq.gz > fastq_cat/10k_PBMC_3p_nextgem_Chromium_X_S1_L001_R2_001.fastq.gz
```

### Run cellranger

```
singularity pull docker://docker.io/nfcore/cellranger:7.1.0
cd $INPUT_DIR && singularity exec -B $(pwd):$(pwd) -B /cfs:/cfs docker://docker.io/nfcore/cellranger:7.1.0 cellranger count --localcores=8 --id=10k_PBMC_3p_nextgem_Chromium_X --fastqs=. --transcriptome=/cfs/sturmgre/tmp/cellranger-ref/refdata-gex-GRCh38-2020-A --localmem=64
cd $MFC && singularity exec -B $(pwd):$(pwd) -B /cfs:/cfs docker://docker.io/nfcore/cellranger:7.1.0 cellranger count --localcores=8 --id=10k_PBMC_3p_nextgem_Chromium_X --fastqs=fc1,fc2,fc3,fc4 --transcriptome=/cfs/sturmgre/tmp/cellranger-ref/refdata-gex-GRCh38-2020-A --localmem=64
cd $MPX && singularity exec -B $(pwd):$(pwd) -B /cfs:/cfs docker://docker.io/nfcore/cellranger:7.1.0 cellranger count --localcores=8 --id=10k_PBMC_3p_nextgem_Chromium_X --fastqs=.  --transcriptome=/cfs/sturmgre/tmp/cellranger-ref/refdata-gex-GRCh38-2020-A --localmem=64
cd fastq_cat && singularity exec -B $(pwd):$(pwd) -B /cfs:/cfs docker://docker.io/nfcore/cellranger:7.1.0 cellranger count --localcores=8 --id=10k_PBMC_3p_nextgem_Chromium_X --fastqs=.  --transcriptome=/cfs/sturmgre/tmp/cellranger-ref/refdata-gex-GRCh38-2020-A --localmem=64
```

### Check results

```
> md5sum **/outs/*.h5 | sort
2a7a5a022f01cb8d95965980f0d95ca5  10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X/outs/raw_feature_bc_matrix.h5
2a7a5a022f01cb8d95965980f0d95ca5  fastq_cat/10k_PBMC_3p_nextgem_Chromium_X/outs/raw_feature_bc_matrix.h5
2a7a5a022f01cb8d95965980f0d95ca5  fastq_multiflowcell/10k_PBMC_3p_nextgem_Chromium_X/outs/raw_feature_bc_matrix.h5
2a7a5a022f01cb8d95965980f0d95ca5  fastq_multiplexed/10k_PBMC_3p_nextgem_Chromium_X/outs/raw_feature_bc_matrix.h5
533f4156f2d90152594ebc587b7fde0a  10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X/outs/filtered_feature_bc_matrix.h5
533f4156f2d90152594ebc587b7fde0a  fastq_cat/10k_PBMC_3p_nextgem_Chromium_X/outs/filtered_feature_bc_matrix.h5
533f4156f2d90152594ebc587b7fde0a  fastq_multiflowcell/10k_PBMC_3p_nextgem_Chromium_X/outs/filtered_feature_bc_matrix.h5
533f4156f2d90152594ebc587b7fde0a  fastq_multiplexed/10k_PBMC_3p_nextgem_Chromium_X/outs/filtered_feature_bc_matrix.h5
cd33d3aa95c0d5cda7cf50908c5399c1  10k_PBMC_3p_nextgem_Chromium_X_fastqs/10k_PBMC_3p_nextgem_Chromium_X/outs/molecule_info.h5
cd33d3aa95c0d5cda7cf50908c5399c1  fastq_cat/10k_PBMC_3p_nextgem_Chromium_X/outs/molecule_info.h5
cd33d3aa95c0d5cda7cf50908c5399c1  fastq_multiflowcell/10k_PBMC_3p_nextgem_Chromium_X/outs/molecule_info.h5
cd33d3aa95c0d5cda7cf50908c5399c1  fastq_multiplexed/10k_PBMC_3p_nextgem_Chromium_X/outs/molecule_info.h5
```
