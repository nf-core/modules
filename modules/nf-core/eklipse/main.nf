
process EKLIPSE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::eklipse=1.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eklipse:1.8--hdfd78af_1':
        'biocontainers/eklipse:1.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam)
    path ref_gb

    output:
    tuple val(meta), path("*deletions.csv"), path("*genes.csv"), emit: csv
    tuple val(meta), path("*.png")                             , emit: circos
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_gb = ref_gb ? "$ref_gb" : "/usr/local/bin/data/NC_012920.1.gb"
    def EKLIPSE_VERSION = "1.8" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "$bam\t${prefix}" > infile.txt
    eKLIPse.py \\
        -in infile.txt \\
        -ref $ref_gb
    mv eKLIPse_*/eKLIPse_deletions.csv eKLIPse_deletions.csv
    mv eKLIPse_*/eKLIPse_genes.csv eKLIPse_genes.csv
    mv eKLIPse_*/eKLIPse_${prefix}.png eKLIPse_${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eKLIPse: $EKLIPSE_VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def EKLIPSE_VERSION = "1.8"
    """
    touch eKLIPse_deletions.csv
    touch eKLIPse_genes.csv
    touch eKLIPse_${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eKLIPse: $EKLIPSE_VERSION
    END_VERSIONS
    """
}
