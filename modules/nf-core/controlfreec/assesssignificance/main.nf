process CONTROLFREEC_ASSESSSIGNIFICANCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d4f13d580f9baf1033ec932a671a3ffaa9df2ed423465d97fab604ecc12b5358/data'
:         'community.wave.seqera.io/library/control-freec:11.6--44126e4b3b49c589' }"

    input:
    tuple val(meta), path(cnvs), path(ratio)

    output:
    tuple val(meta), path("*.p.value.txt"), emit: p_value_txt
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('controlfreec'), val("11.6b"), emit: versions_controlfreec, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat \$(which assess_significance.R) | R --slave --args ${cnvs} ${ratio}

    mv *.p.value.txt ${prefix}.p.value.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.p.value.txt
    """
}
