process LJA {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/lja:${task.ext.version ?: '0.2--h5b5514e_2'}": "biocontainers/lja:${task.ext.version ?: '0.2--h5b5514e_2'}" }"
    container "docker://troder/lja"

    input:
    tuple val(meta), path(reads)
//    val ploidy

    output:
    // tuple val(meta), path("*.fasta"),    emit: fastaraw
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.gfa.gz")  , emit: gfa
    tuple val(meta), path("*.stdout")  , emit: stdout
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    def input_list = reads.collect{"--reads $it"}.join(' ')

    """
    lja \\
       $args \\
        $input_list \\
        --output-dir . \\
        --threads $task.cpus \\
        > ${prefix}.lja.stdout

    gzip -c -n assembly.fasta > ${prefix}.fasta.gz
    gzip -c -n mdbg.gfa > ${prefix}.gfa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LJA: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    echo stub | gzip -c > ${prefix}.fasta.gz
    echo stub | gzip -c > ${prefix}.gfa.gz

    cat <<-END_STDOUT > ${prefix}.lja.stdout
    00:00:00 0Mb  INFO: LJA pipeline finished
    END_STDOUT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LJA: $VERSION
    END_VERSIONS
    """
}
