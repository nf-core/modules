#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.outdir = "."
params.genome = ""
params.bowtie2_args = ''
// Bowtie2 arguments should be supplied in the following format to work:
// --bowtie2_args="--score-min L,0,-0.8"

params.verbose = false

if (params.verbose){
    println ("[WORKFLOW] BOWTIE2 ARGS: "      + params.bowtie2_args)
}

// for other genomes this needs to be handled somehow to return all possible genomes
genomeValues = ["name" : params.genome]
genomeValues["bowtie2"] = "/bi/home/fkrueger/VersionControl/nf-core-modules/test-datasets/indices/bowtie2/E_coli/${params.genome}";

include '../main.nf'   params(genome: genomeValues)

ch_read_files = Channel
  .fromFilePairs('../../../test-datasets/Ecoli*{1,2}.fastq.gz',size:-1)
  // .view()  // to check whether the input channel works

workflow {

    main:
        BOWTIE2(ch_read_files, params.outdir, params.bowtie2_args, params.verbose)
    
}