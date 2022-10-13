process ULTRA_INDEX {
    tag "$gtf"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ultra_bioinformatics=0.0.4" : null)
    container { workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.0.4.1--pyh5e36f6f_0':
        "${params.docker_url ?: 'quay.io/biocontainers'}/ultra_bioinformatics:0.0.4.1--pyh5e36f6f_0" }

    input:
    path fasta
    path gtf

    output:
    tuple path("*.pickle"), path("*.db"), emit: index
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    uLTRA \\
        index \\
        $args \\
        $fasta \\
        $gtf \\
        ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultra: \$( uLTRA --version|sed 's/uLTRA //g' )
    END_VERSIONS
    """
}
