#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.outdir = "."
params.genome = ""
params.hisat2_args = ''
// HISAT2 arguments should be supplied in the following format to work:
// --hisat2_args="--score-min L,0,-0.8"

params.verbose = false

if (params.verbose){
    println ("[WORKFLOW] HISAT2 ARGS ARE: "       + params.hisat2_args)
}
// for other genomes this needs to be handled somehow to return all possible genomes
genomeValues = ["name" : params.genome]
genomeValues["hisat2"] = "/bi/home/fkrueger/VersionControl/nf-core-modules/test-datasets/indices/hisat2/E_coli/${params.genome}";

include '../main.nf'   params(genome: genomeValues)

ch_read_files = Channel
  .fromFilePairs('../../../test-datasets/Ecoli*{1,2}.fastq.gz',size:-1)
  // .view()  // to check whether the input channel works

workflow {

    main:
        HISAT2(ch_read_files, params.outdir, params.hisat2_args, params.verbose)
}





