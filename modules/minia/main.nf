process MINIA {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::minia=3.2.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minia:3.2.6--h9a82719_0' :
        'quay.io/biocontainers/minia:3.2.6--h9a82719_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.contigs.fa'), emit: contigs
    tuple val(meta), path('*.unitigs.fa'), emit: unitigs
    tuple val(meta), path('*.h5')        , emit: h5
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_list = reads.join(",")
    """
    echo "${read_list}" | sed 's/,/\\n/g' > input_files.txt
    minia \\
        $args \\
        -nb-cores $task.cpus \\
        -in input_files.txt \\
        -out $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minia: \$(echo \$(minia --version 2>&1 | grep Minia) | sed 's/^.*Minia version //;')
    END_VERSIONS
    """
}
