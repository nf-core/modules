process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(intervals), val(scale)
    path  sizes
    val   extension

    output:
    tuple val(meta), path("*.${extension}"), emit: genomecov
    path  "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args_list = args.tokenize()
    args += (scale > 0 && scale != 1) ? " -scale $scale" : ""
    if (!args_list.contains('-bg') && (scale > 0 && scale != 1)) {
        args += " -bg"
    }

    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (intervals.name =~ /\.bam/) {
        """
        bedtools \\
            genomecov \\
            -ibam $intervals \\
            $args \\
            > ${prefix}.${extension}

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        bedtools \\
            genomecov \\
            -i $intervals \\
            -g $sizes \\
            $args \\
            > ${prefix}.${extension}

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
}
