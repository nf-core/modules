process RIBOCODE_RIBOCODE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocode:1.2.15--pyhfa5458b_0':
        'biocontainers/ribocode:1.2.15--pyhfa5458b_0' }"

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
    # Run RiboCode and capture output using tee to both display and save
    RiboCode \\
        -a $annotation \\
        -c $config \\
        -o ${prefix} \\
        $args 2>&1 | tee ribocode_output.log || true

    # Check if RiboCode output contains any error messages
    if grep -qi "^Error" ribocode_output.log; then
        echo "" >&2
        echo "ERROR: RiboCode failed. Check the output above for details." >&2
        echo "Common causes:" >&2
        echo "  - Invalid config file from metaplots (try adjusting --extra_ribocode_metaplots_args '-f0_percent 0.XX')" >&2
        echo "  - Insufficient data from Ribo-Seq alignment" >&2
        exit 1
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt
    touch ${prefix}_collapsed.txt
    touch ${prefix}_ORFs_category.pdf
    touch ${prefix}_psites.hd5
    """
}
