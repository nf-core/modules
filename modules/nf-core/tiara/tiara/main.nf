process TIARA_TIARA {
    tag "${meta.id}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5909e63835fda5723b1408001568a306e234dd0035b3518af8a8930625e449d/data'
        : 'community.wave.seqera.io/library/tiara:1.0.3--e367cb0eaef39fa3'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.{txt,txt.gz}"), emit: classifications
    tuple val(meta), path("log_*.{txt,txt.gz}")    , emit: log
    tuple val(meta), path("*.{fasta,fasta.gz}")    , emit: fasta, optional: true
    tuple val("${task.process}"), val('tiara'), val("1.0.3"), topic: versions, emit: versions_tiara
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    tiara -i ${fasta} \
        -o ${prefix}.txt \
        --threads ${task.cpus} \
        ${args}

    ## fix gzip flag weirdness and ensure consistent .fasta filename output
    ## check if fasta files are being output
    if echo "${args}" | grep -qE "tf|to-fasta"; then
        ## check if we've asked for gzip output, then rename files consistently
        if echo "${args}" | grep -q "gz"; then
            find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}.fasta.gz"' \\;
        else
            find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}.fasta"' \\;
        fi
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch log_${prefix}.txt
    touch bacteria_${prefix}.fasta
    """
}
