process FLYE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        'biocontainers/flye:2.9--py39h6935b12_1' }"

    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.gfa.gz")  , emit: gfa
    tuple val(meta), path("*.gv.gz")   , emit: gv
    tuple val(meta), path("*.txt")     , emit: txt
    tuple val(meta), path("*.log")     , emit: log
    tuple val(meta), path("*.json")    , emit: json
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Flye. Options: ${valid_mode.join(', ')}" }
    """
    flye \\
        $mode \\
        $reads \\
        --out-dir . \\
        --threads \\
        $task.cpus \\
        $args

    gzip -c assembly.fasta > ${prefix}.assembly.fasta.gz
    gzip -c assembly_graph.gfa > ${prefix}.assembly_graph.gfa.gz
    gzip -c assembly_graph.gv > ${prefix}.assembly_graph.gv.gz
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv flye.log ${prefix}.flye.log
    mv params.json ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub | gzip -c > ${prefix}.assembly.fasta.gz
    echo stub | gzip -c > ${prefix}.assembly_graph.gfa.gz
    echo stub | gzip -c > ${prefix}.assembly_graph.gv.gz
    echo contig_1 > ${prefix}.assembly_info.txt
    echo stub > ${prefix}.flye.log
    echo stub > ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """
}
