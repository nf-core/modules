process CELESTA {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/schapirolabor/mcmicro-celesta:v0.0.2"

    input:
    tuple val(meta), path(img_data)
    path(signature)
    path(high_thresholds)
    path(low_thresholds)

    output:
    tuple val(meta), path("*results.csv"), emit: celltypes
    path "*quality.csv"                  , emit: quality
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "celesta module was created only for Docker, Singularity or Podman. It does not support Conda!"
    }
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def low_thresholds_cmd = low_thresholds ? "--low $low_thresholds" : ""
    def VERSION = '1.0.0'

    """
    Rscript /local/CELESTA_CLI.R \\
        -i $img_data \\
        -s $signature \\
        --high $high_thresholds \\
        $low_thresholds_cmd \\
        -o . \\
        -t $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        celesta: $VERSION
    END_VERSIONS
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "celesta module was created only for Docker, Singularity or Podman. It does not support Conda!"
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'

    """
    touch ${prefix}_celesta_stub_results.csv
    touch ${prefix}_stub_quality.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        celesta: $VERSION
    END_VERSIONS
    """
}
