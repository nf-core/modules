nextflow.preview.dsl = 2

process FASTQC {

    // tag "FastQC - $sample_id"

    input:
        tuple val(name), path(reads)
        val (outputdir)
        // fastqc_args are best passed into the workflow in the following manner:
        // --fastqc_args="--nogroup -a custom_adapter_file.txt"
        val (fastqc_args)
        val (verbose)

    output:
        tuple val(name), path ("*fastqc*"), emit: all
        path "*.zip",                       emit: report // e.g. for MultiQC later

    // container 'quay.io/biocontainers/fastqc:0.11.8--2'

    publishDir "$outputdir",
        mode: "copy", overwrite: true

    script:

        if (verbose){
            println ("[MODULE] FASTQC ARGS: " + fastqc_args)
        }

        """
        module load fastqc
        fastqc $fastqc_args -q -t 2 $reads

        fastqc --version &> fastqc.version.txt
        """

}
