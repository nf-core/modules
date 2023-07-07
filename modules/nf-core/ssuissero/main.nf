process SSUISSERO {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::ssuissero=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ssuissero%3A1.0.1--hdfd78af_0':
        'biocontainers/ssuissero:1.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    SsuisSero.sh \\
        -i $fasta_name \\
        -o ./ \\
        -s $prefix \\
        -x fasta \\
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ssuissero: $VERSION
    END_VERSIONS
    """
}
