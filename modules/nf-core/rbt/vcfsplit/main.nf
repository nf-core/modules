process RBT_VCFSPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rust-bio-tools:0.42.2--h4458251_1':
        'biocontainers/rust-bio-tools:0.42.2--h4458251_1' }"

    input:
    tuple val(meta), path(vcf)
    val(numchunks)

    output:
    tuple val(meta), path("*"), emit: bcfchunks
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunks = task.ext.numchunks ? (task.ext.numchunks - 1) : 15
    """
    rbt vcf-split ${vcf} ${prefix}{0..3}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rust-Bio-Tools: \$$(rbt --version | grep -o '[0-9.]\+')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunks = task.ext.numchunks ? (task.ext.numchunks - 1) : 15
    """
    touch hell.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rust-Bio-Tools: \$(rbt --version | grep -o '[0-9.]\+')
    END_VERSIONS
    """
}
