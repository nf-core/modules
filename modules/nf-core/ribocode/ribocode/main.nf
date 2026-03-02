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
    set -o pipefail

    RiboCode \\
        -a $annotation \\
        -c $config \\
        -o ${prefix} \\
        $args 2>&1 | tee ${prefix}_ribocode.log

    # RiboCode may exit 0 even on failure (e.g. P-site detection errors).
    # Verify expected output was actually produced.
    if [ ! -s ${prefix}.txt ]; then
        echo "ERROR: RiboCode produced no output - check log for details." >&2
        echo "A common cause is failed P-site detection in metaplots." >&2
        echo "Try relaxing metaplots thresholds via extra_ribocode_metaplots_args" >&2
        echo "  e.g. '-pv1 0.05 -pv2 0.05 -f0_percent 0.5'" >&2
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
