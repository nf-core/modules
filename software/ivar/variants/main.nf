// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/1.3.1--h089eab3_0"
    } else {
        container "quay.io/biocontainers/1.3.1--h089eab3_0"
    }

    input:
    tuple val(meta), path(bam)
    path(reference_fasta)
    path(gff_file)

    output:
    tuple val(meta), path("*.tsv"), emit: variants
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // the gff file is optional, so following the pattern suggested here:
    //https://github.com/nextflow-io/patterns/blob/master/docs/optional-input.adoc
    def gff      = gff_file.name != 'NO_FILE' ? "-g $gff_file" : ""
    """
    samtools mpileup \\
        -aa -A -d 0 -B -Q 0 \\
        $options.args2 \\
        --reference $reference_fasta \\
        $bam | \\
        ivar variants \\
        -p $prefix \\
        $options.args \\
        -r $reference_fasta \\
        $gff
    """
}
