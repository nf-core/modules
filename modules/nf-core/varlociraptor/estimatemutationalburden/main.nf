process VARLOCIRAPTOR_ESTIMATEMUTATIONALBURDEN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c99922d373e15a47f7a4ad4692e859f99e97a00c7872214903e655f66019c772/data'
        : 'community.wave.seqera.io/library/varlociraptor:8.8.0--0d1bd46f0f5e39af'}"

    input:
    tuple val(meta), path(vcf)
    val(output_mode)

    output:
    tuple val(meta), path("*.tsv"), emit: table, optional: true
    tuple val(meta), path("*.svg"), emit: svg,   optional: true
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = output_mode == 'curve' ? "| vg2svg > ${prefix}.svg" : "> ${prefix}.tsv"
    """
    varlociraptor estimate mutational-burden --coding-genome-size 3e7 --events SOMATIC_TUMOR --sample tumor < calls.bcf > tmb.tsv
    varlociraptor estimate mutational-burden \\
        --coding-genome-size 3e7 --events SOMATIC_TUMOR --sample tumor \\
        --mode ${output_mode} \\
        ${args} \\
        < ${vcf} \\
        ${output_cmd}

    // --coding-genome-size 3e7 --events SOMATIC_TUMOR --sample tumor
    // probably needs to be added

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """
}
