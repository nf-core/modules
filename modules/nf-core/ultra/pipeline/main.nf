process ULTRA_PIPELINE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ultra_bioinformatics=0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.1--pyh7cba7a3_1':
        'biocontainers/ultra_bioinformatics:0.1--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(reads)
    path genome
    path gtf

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    uLTRA \\
        pipeline \\
        --t $task.cpus \\
        --prefix $prefix \\
        $args \\
        $genome \\
        $gtf \\
        $reads \\
        ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
    END_VERSIONS
    """
}
