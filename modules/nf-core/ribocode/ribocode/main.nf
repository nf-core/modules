process RIBOCODE_RIBOCODE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe815db0864b45b91afc7bc84c55cb60acb0035e7248dda7f480a55c4cb105d7/data':
        'community.wave.seqera.io/library/ribocode:1.2.15--5530b252f5433a62' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(annotation)
    tuple val(meta3), path(config)

    output:

    tuple val(meta), path("*.txt")                                                  , emit: orf_txt
    tuple val(meta), path("*_collapsed.txt")                                        , emit: orf_txt_collapsed
    tuple val(meta), path("*_ORFs_category.pdf")                                    , emit: orf_pdf, optional: true
    tuple val(meta), path("*_psites.hd5")                                           , emit: psites_hd5, optional: true
    tuple val("${task.process}"), val('ribocode'), eval('RiboCode --version  2>&1') , emit: versions_ribocode, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Run RiboCode and capture output to check for errors
    RiboCode \\
        -a $annotation \\
        -c $config \\
        -o ${prefix} \\
        $args 2>&1 | tee ribocode_output.log || true

    # Check if RiboCode output contains error messages
    if grep -qiE "^Error|Error:" ribocode_output.log; then
        echo "ERROR: RiboCode failed. Check the output above for details." >&2
        echo "Common causes:" >&2
        echo "  - Invalid config file from metaplots (try lowering cutoff via ext.args, e.g., '-f0_percent 0.5')" >&2
        echo "  - Insufficient data from Ribo-Seq alignment" >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt
    touch ${prefix}_collapsed.txt
    touch ${prefix}_ORFs_category.pdf
    touch ${prefix}_psites.hd5
    """
}
