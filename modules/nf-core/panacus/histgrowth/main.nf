process PANACUS_HISTGROWTH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::panacus=0.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/panacus:0.2.3--h031d066_0':
        'biocontainers/panacus:0.2.3--h031d066_0' }"

    input:
    tuple val(meta), path(gfa)
    path(bed_subset)
    path(bed_exclude)
    path(tsv_groupby)

    output:
    tuple val(meta), path("*.{tsv, html}"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-format table") || args.contains("-o table") ? "tsv" :
                    args.contains("--output-format html") || args.contains("-o html") ? "html" :
                    "tsv"
    def subset_query  = bed_subset ? "--subset ${bed_subset}" : ""
    def exclude_query = bed_exclude ? "--exclude ${bed_exclude}" : ""
    def groupby_query = tsv_groupby ? "--groupby ${tsv_groupby}" : ""
    """
    panacus \\
        histgrowth \\
        $args \\
        $subset_query \\
        $exclude_query \\
        $groupby_query \\
        --threads $task.cpus \\
        $gfa > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panacus: \$(echo \$(panacus --version) | sed 's/^panacus //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-format table") || args.contains("-o table") ? "tsv" :
                    args.contains("--output-format html") || args.contains("-o html") ? "html" :
                    "tsv"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panacus: \$(echo \$(panacus --version) | sed 's/^panacus //' ))
    END_VERSIONS
    """
}
