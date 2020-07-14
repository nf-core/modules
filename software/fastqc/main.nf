nextflow.preview.dsl = 2

process FASTQC {
    input:
        tuple val(name), path(reads)
        val (fastqc_args)
        val (outputdir)

    output:
        tuple val(name), path ("*fastqc*"), emit: all
        path "*.zip",                       emit: report // e.g. for MultiQC later

    container 'quay.io/biocontainers/fastqc:0.11.8--2'

    publishDir "$outputdir",
        mode: "copy", overwrite: true

    script:
        """
        fastqc $fastqc_args -q -t 2 $reads
        fastqc --version &> fastqc.version.txt
        """
}
