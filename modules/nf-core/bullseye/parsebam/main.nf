process BULLSEYE_PARSEBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c267028d70907ccfc19fd125a704f7077197254083d97f5aaedda1bdb751c7b/data'
        : 'community.wave.seqera.io/library/bullseye:1.0.0--bf11337413b3307b'}"

    input:
    tuple val(meta), path(bam), path(bai)
    val(bullseye_parse_min_coverage)

    output:
    tuple val(meta), path("*.matrix"), emit: matrix
    tuple val("${task.process}"), val('bullseye'), eval("grep -m1 -i 'version' \$(command -v parseBAM.pl) | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?'"), topic: versions, emit: versions_bullseye

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    perl -S parseBAM.pl \\
        --input ${bam} \\
        --output ${prefix}.matrix \\
        --cpu ${task.cpus} \\
        --minCoverage ${bullseye_parse_min_coverage} \\
        $args

    # parseBAM.pl may emit either <prefix>.matrix or <prefix>.matrix.gz depending on mode/version.
    # Compression and tabix-indexing are left to a dedicated downstream step (e.g. htslib/bgziptabix).
    if [[ -s ${prefix}.matrix.gz ]]; then
        gunzip -f ${prefix}.matrix.gz
    elif [[ ! -s ${prefix}.matrix ]]; then
        echo "parseBAM.pl did not produce ${prefix}.matrix or ${prefix}.matrix.gz" >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.matrix
    """
}
