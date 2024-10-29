
process EKLIPSE {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eklipse:1.8--hdfd78af_1':
        'biocontainers/eklipse:1.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_gb

    output:
    tuple val(meta), path("*deletions.csv") , emit: deletions
    tuple val(meta), path("*genes.csv")     , emit: genes
    tuple val(meta), path("*.png")          , emit: circos
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_gb = ref_gb ? "$ref_gb" : "/usr/local/bin/data/NC_012920.1.gb"
    def VERSION = "1.8" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "$bam\t${prefix}" > infile.txt
    eKLIPse.py \\
        -in infile.txt \\
        $args \\
        -ref $ref_gb
    mv eKLIPse_*/eKLIPse_deletions.csv eKLIPse_${prefix}_deletions.csv
    mv eKLIPse_*/eKLIPse_genes.csv eKLIPse_${prefix}_genes.csv
    mv eKLIPse_*/eKLIPse_${prefix}.png eKLIPse_${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eklipse: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.8"
    """
    touch eKLIPse_${prefix}_deletions.csv
    touch eKLIPse_${prefix}_genes.csv
    touch eKLIPse_${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eklipse: $VERSION
    END_VERSIONS
    """
}
