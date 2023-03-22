process CUTESV {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::cutesv=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutesv:1.0.12--pyhdfd78af_0' :
        'quay.io/biocontainers/cutesv:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)

    output:
    tuple val(meta), path("*_cuteSV.vcf"), emit: sv_vcf // vcf files
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genotyping = params.enable_genotyping ? "--genotyping" : ''
    """
    cuteSV \
        ${bam} \\
        ${fasta} \\
        ${meta.id}_cuteSV.vcf \\
        . \\
        --threads $task.cpus \\
        --sample ${meta.id} \\
        $genotyping \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuteSV: \$( cuteSV --version 2>&1 | sed 's/cuteSV //g' )
    END_VERSIONS
    """
}
