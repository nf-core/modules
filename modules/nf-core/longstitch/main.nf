process LONGSTITCH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::longstitch=1.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longstitch:1.0.5--hdfd78af_0':
        'biocontainers/longstitch:1.0.5--hdfd78af_0' }"

    input:
    tuple val(meta) , path (assembly, stageAs: "assembly.fa")
    tuple val(meta2), path (reads   , stageAs: "reads.fq.gz")
    val mode
    val genome_size

    output:
    tuple val(meta), path("*.fa"), emit: assembly
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def valid_mode = ["tigmint-ntLink-arks", "tigmint-ntLink", "ntLink-arks"]
    def genome_size_value = genome_size ? "G=$genome_size" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    longstitch \\
        $mode \\
        $args \\
        $genome_size_value \\
        t=$task.cpus \\
        draft=assembly \\
        reads=reads \\
        out_prefix=${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longstitch : \$( longstitch  --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longstitch : \$( longstitch  --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
