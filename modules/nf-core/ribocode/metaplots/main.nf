process RIBOCODE_METAPLOTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe815db0864b45b91afc7bc84c55cb60acb0035e7248dda7f480a55c4cb105d7/data':
        'community.wave.seqera.io/library/ribocode:1.2.15--5530b252f5433a62' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(annotation)

    output:
    tuple val(meta), path("*config.txt")                                            , emit: config
    tuple val(meta), path("*.pdf")                                                  , emit: pdf
    tuple val("${task.process}"), val('ribocode'), eval('RiboCode --version  2>&1') , emit: versions_ribocode, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metaplots \\
        -a $annotation \\
        -r $bam \\
        -o ${prefix} \\
        $args

    # Check config file has a sample config line (non-empty, doesn't start with #)
    if ! grep -qE '^[^#[:space:]]' ${prefix}_pre_config.txt; then
        echo "ERROR: metaplots created config file with no data (only header)." >&2
        echo "This usually indicates insufficient periodic signal in Ribo-Seq data." >&2
        echo "Consider lowering the cutoff via ext.args (e.g., '-f0_percent 0.5')." >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_config.txt
    touch ${prefix}_report.pdf
    """
}
