process BULLSEYE_PARSEBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d1b3ef00779596d2255fc5dffeae072301ccc37d6cdb941c2aa8d67534ced41/data'
        : 'community.wave.seqera.io/library/bullseye:1.1.0--6a0a52edc73ac840'}"

    input:
    tuple val(meta), path(bam), path(bai)
    val(bullseye_parse_min_coverage)

    output:
    tuple val(meta), path("*.matrix.gz"), emit: matrix
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
    # Tabix-indexing is left to a dedicated downstream step (e.g. htslib/bgziptabix).
    if [[ -s ${prefix}.matrix ]]; then
        bgzip -f ${prefix}.matrix
    elif [[ ! -s ${prefix}.matrix.gz ]]; then
        echo "parseBAM.pl did not produce ${prefix}.matrix or ${prefix}.matrix.gz" >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.matrix.gz
    """
}
