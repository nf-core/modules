process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(intervals), val(scale)
    path  sizes
    val   extension

    output:
    tuple val(meta), path("*.${extension}"), emit: genomecov
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_list = args.tokenize()
    args += (scale > 0 && scale != 1) ? " -scale $scale" : ""
    if (!args_list.contains('-bg') && (scale > 0 && scale != 1)) {
        args += " -bg"
    }

    def prefix = task.ext.prefix ?: "${meta.id}"
    if (intervals.name =~ /\.bam/) {
        """
        bedtools \\
            genomecov \\
            -ibam $intervals \\
            $args \\
            > ${prefix}.${extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
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
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
