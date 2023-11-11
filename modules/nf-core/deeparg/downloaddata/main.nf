process DEEPARG_DOWNLOADDATA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeparg:1.0.2--pyhdfd78af_1' :
        'biocontainers/deeparg:1.0.2--pyhdfd78af_1' }"

    input:

    output:
    path "db/"          , emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def DATA_URL='https://zenodo.org/records/8280582/files/deeparg.zip' // As per https://github.com/gaarangoa/deeparg/blob/537f0394daf4de85858a390e82a3d833febe280d/deeparg/entry.py#L77
    def VERSION='1.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    wget -O deeparg.zip $DATA_URL
    unzip deeparg.zip
    mv deeparg db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeparg: $VERSION
    END_VERSIONS
    """
}
