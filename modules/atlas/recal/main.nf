process ATLAS_RECAL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::atlas=0.9.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/atlas:0.9.9--h082e891_0':
        'quay.io/biocontainers/atlas:0.9.9--h082e891_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(empiric)
    tuple path(fasta), path(alleles)

    output:
    tuple val(meta), path("*.txt"), emit:recal_patterns
    path "versions.yml"           , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def PMD = empiric ? "pmdFile=${empiric}" : ""

    """
    atlas \\
        task=recal \\
        bam=$bam \\
        $PMD \\
        regions=$alleles \\
        out=$prefix \\
        $args 


    (atlas 2>&1) | grep Atlas | head -n 1 | sed -e 's/^[ \t]*Atlas //' > versions.yml


    """
}
