process CUSTOM_ORFNORMALISE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15f690593e9d2bd7ebeb72fff4068eeae9cdfc909103de457b81a25d29cfbef3/data' :
        'community.wave.seqera.io/library/python_pyyaml:25b0f22c7e7bf847' }"

    input:
    tuple val(meta) , path(orfs_table), val(caller)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("${prefix}.bed12"), emit: bed12
    tuple val(meta), path("${prefix}.tsv")  , emit: tsv
    path "versions.yml"                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix ?: "${meta.id}.normalised"
    sample_id = meta.id ?: 'unknown'
    args      = task.ext.args ?: ''
    template 'orfnormalise.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.normalised"
    """
    touch ${prefix}.bed12
    touch ${prefix}.tsv

    python -c "import platform, yaml; yaml.safe_dump({'${task.process}': {'python': platform.python_version()}}, open('versions.yml', 'w'), default_flow_style=False, sort_keys=False)"
    """
}
