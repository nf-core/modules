process INTEGRONFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/integron_finder:2.0.5--pyhdfd78af_0':
        'biocontainers/integron_finder:2.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*/*.gbk")      , emit: gbk, optional: true
    tuple val(meta), path("*/*.integrons"), emit: integrons
    tuple val(meta), path("*/*.summary")  , emit: summary
    path("*/integron_finder.out")         , emit: out
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args ?: ''
    def is_compressed_fasta = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name          = fasta.getName().replace(".gz", "")
    
    """
    if [ "$is_compressed_fasta" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    integron_finder \\
        $args \\
        --local-max \\
        --cpu $task.cpus \\
        --promoter-attI \\
        --eagle-eyes \\
        --gbk \\
        $fasta_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integronfinder: \$(integron_finder --version | sed -n 's/.*version \\(.*\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p Results_Integron_Finder_${prefix}
    touch "Results_Integron_Finder_${prefix}/${prefix}.gbk"
    touch "Results_Integron_Finder_${prefix}/${prefix}.integrons"
    touch "Results_Integron_Finder_${prefix}/${prefix}.summary"
    touch "Results_Integron_Finder_${prefix}/integron_finder.out"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integronfinder: \$(integron_finder --version | sed -n 's/.*version \\(.*\\).*/\\1/p')
    END_VERSIONS
    """
}