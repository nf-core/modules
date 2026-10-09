process KRAKEN2_BUILDSTANDARD {
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/ec5af88da27c32c52b667af2af67ce76ffc96715c1a40e0f4ab72dd87623adfa/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:b111fc3860b64b23' }"

    input:
    val cleaning

    output:
    path ("${prefix}"), emit: db
    tuple val("${task.process}"), val('kraken2'), eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"'), topic: versions, emit: versions_kraken2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "kraken2_standard_db"
    runclean = cleaning ? "kraken2-build --clean --db ${prefix}" : ""
    """
    kraken2-build \\
        --standard \\
        ${args} \\
        --threads ${task.cpus} \\
        --db ${prefix}

    ${runclean}
    """

    stub:
    prefix = task.ext.prefix ?: "kraken2_standard_db"
    """
    mkdir -p "${prefix}"
    """
}
