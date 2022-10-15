process ULTRA_PIPELINE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::ultra_bioinformatics=0.0.4.1" : null)
    def container_image = "ultra_bioinformatics:0.0.4.1--pyh5e36f6f_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


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
