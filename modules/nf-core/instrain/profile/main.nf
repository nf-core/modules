process INSTRAIN_PROFILE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::instrain=1.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/instrain:1.6.1--pyhdfd78af_0':
        'biocontainers/instrain:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path genome_fasta
    path genes_fasta
    path stb_file

    output:
    tuple val(meta), path("*.IS") , emit: profile
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes_args = genes_fasta ? "-g ${genes_fasta}": ''
    def stb_args = stb_file ? "-s ${stb_file}": ''
    """
    inStrain \\
        profile \\
        $bam \\
        $genome_fasta \\
        -o ${prefix}.IS \\
        -p $task.cpus \\
        $genes_args \\
        $stb_args \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(echo \$(inStrain profile --version 2>&1) | awk 'NF{ print \$NF }')
    END_VERSIONS
    """
}
