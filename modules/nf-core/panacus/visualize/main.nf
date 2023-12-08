process PANACUS_VISUALIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/panacus:0.2.3--h031d066_0':
        'biocontainers/panacus:0.2.3--h031d066_0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.{eps,jpg,jpeg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff,webp}"), emit: image
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--format eps") || args.contains("-f eps") ? "eps" :
                    args.contains("--format jpg") || args.contains("-f jpg") ? "jpg" :
                    args.contains("--format jpeg") || args.contains("-f jpeg") ? "jpeg" :
                    args.contains("--format pdf") || args.contains("-f pdf") ? "pdf" :
                    args.contains("--format pgf") || args.contains("-f pgf") ? "pgf" :
                    args.contains("--format png") || args.contains("-f png") ? "png" :
                    args.contains("--format ps") || args.contains("-f ps") ? "ps" :
                    args.contains("--format raw") || args.contains("-f raw") ? "raw" :
                    args.contains("--format rgba") || args.contains("-f rgba") ? "rgba" :
                    args.contains("--format svg") || args.contains("-f svg") ? "svg" :
                    args.contains("--format svgz") || args.contains("-f svgz") ? "svgz" :
                    args.contains("--format tif") || args.contains("-f tif") ? "tif" :
                    args.contains("--format tiff") || args.contains("-f tiff") ? "tiff" :
                    args.contains("--format webp") || args.contains("-f webp") ? "webp" :
                    "pdf"
    def output_pipe = args.contains("--split_subfigures") ? "" : "> ${prefix}.${extension}"
    """
    panacus-visualize \\
        $args \\
        $tsv \\
        $output_pipe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panacus: \$(echo \$(panacus --version) | sed 's/^panacus //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--format eps") || args.contains("-f eps") ? "eps" :
                    args.contains("--format jpg") || args.contains("-f jpg") ? "jpg" :
                    args.contains("--format jpeg") || args.contains("-f jpeg") ? "jpeg" :
                    args.contains("--format pdf") || args.contains("-f pdf") ? "pdf" :
                    args.contains("--format pgf") || args.contains("-f pgf") ? "pgf" :
                    args.contains("--format png") || args.contains("-f png") ? "png" :
                    args.contains("--format ps") || args.contains("-f ps") ? "ps" :
                    args.contains("--format raw") || args.contains("-f raw") ? "raw" :
                    args.contains("--format rgba") || args.contains("-f rgba") ? "rgba" :
                    args.contains("--format svg") || args.contains("-f svg") ? "svg" :
                    args.contains("--format svgz") || args.contains("-f svgz") ? "svgz" :
                    args.contains("--format tif") || args.contains("-f tif") ? "tif" :
                    args.contains("--format tiff") || args.contains("-f tiff") ? "tiff" :
                    args.contains("--format webp") || args.contains("-f webp") ? "webp" :
                    "pdf"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panacus: \$(echo \$(panacus --version) | sed 's/^panacus //' ))
    END_VERSIONS
    """
}
