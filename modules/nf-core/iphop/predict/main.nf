process IPHOP_PREDICT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::iphop=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iphop:1.3.1--pyhdfd78af_0':
        'biocontainers/iphop:1.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path  iphop_db

    output:
    tuple val(meta), path("Host_prediction_to_genus_m*.csv"), emit: iphop_genus

    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    iphop \\
        predict \\
        --fa_file $fasta \\
        --out_dir . \\
        --db_dir $iphop_db \\
        --num_threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """
}
