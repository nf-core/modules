process ADAPTERREMOVAL {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::adapterremoval=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/adapterremoval:2.3.2--hb7ba0dd_0' :
        'quay.io/biocontainers/adapterremoval:2.3.2--hb7ba0dd_0' }"

    input:
    tuple val(meta), path(reads)
    path(adapterlist)

    output:
    tuple val(meta), path("${prefix}.truncated.gz")            , optional: true, emit: singles_truncated
    tuple val(meta), path("${prefix}.discarded.gz")            , optional: true, emit: discarded
    tuple val(meta), path("${prefix}.pair1.truncated.gz")      , optional: true, emit: pair1_truncated
    tuple val(meta), path("${prefix}.pair2.truncated.gz")      , optional: true, emit: pair2_truncated
    tuple val(meta), path("${prefix}.collapsed.gz")            , optional: true, emit: collapsed
    tuple val(meta), path("${prefix}.collapsed.truncated.gz")  , optional: true, emit: collapsed_truncated
    tuple val(meta), path("${prefix}.paired.gz")               , optional: true, emit: paired_interleaved
    tuple val(meta), path('*.log')                             , emit: log
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def list = adapterlist ? "--adapter-list ${adapterlist}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        AdapterRemoval  \\
            --file1 $reads \\
            $args \\
            $adapterlist \\
            --basename ${prefix} \\
            --threads ${task.cpus} \\
            --settings ${prefix}.log \\
            --seed 42 \\
            --gzip

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            adapterremoval: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    } else {
        """
        AdapterRemoval  \\
            --file1 ${reads[0]} \\
            --file2 ${reads[1]} \\
            $args \\
            $adapterlist \\
            --basename ${prefix} \\
            --threads $task.cpus \\
            --settings ${prefix}.log \\
            --seed 42 \\
            --gzip

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            adapterremoval: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    }

}
