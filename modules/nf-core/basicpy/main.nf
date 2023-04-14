process BASICPY {
    tag "$meta.id"
    label 'process_single'

    container "yfukai/basicpy-docker-mcmicro:0.1.2"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.tiff"), emit: fields
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.1.2"
    """
    /opt/main.py --cpu $image .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basicpy:: $VERSION
    END_VERSIONS
    """
}
