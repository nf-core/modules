process VARLOCIRAPTOR_CALLVARIANTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c99922d373e15a47f7a4ad4692e859f99e97a00c7872214903e655f66019c772/data'
        : 'community.wave.seqera.io/library/varlociraptor:8.8.0--0d1bd46f0f5e39af'}"

    input:
    tuple val(meta), path(vcfs)
    path scenario
    val scenario_aliases

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_called"

    //If we use a scenario file and if there is more than 1 vcf, then collect scenario_aliaes and vcfs to scenario_alias_0=vcf_0 scenario_alias_1=vcf_1, etc
    //If we use a scenario file and if there is exactly 1 vcf, then scenario_alias=vcf
    def scenario_samples = vcfs instanceof List && vcfs.size() > 1 ? [scenario_aliases, vcfs].transpose().collect { "${it[0]}=${it[1]}" }.join(' ') : "${scenario_aliases}=${vcfs}"
    """
    varlociraptor call variants \\
        --output ${prefix}.bcf \\
        generic --scenario ${scenario} --obs ${scenario_samples} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_called"
    """
    touch ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """
}
