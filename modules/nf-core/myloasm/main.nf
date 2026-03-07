process MYLOASM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/myloasm:0.2.0--ha6fb395_0':
        'community.wave.seqera.io/library/myloasm:0.2.0--036e61a36965d08c' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}")                                           , emit: results
    tuple val(meta), path("${prefix}/assembly_primary.fa")                       , emit: contigs
    tuple val(meta), path("${prefix}/final_contig_graph.gfa")                    , emit: gfa
    tuple val(meta), path("${prefix}/alternate_assemblies/assembly_alternate.fa"), emit: contigs_alt
    tuple val(meta), path("${prefix}/alternate_assemblies/duplicated_contigs.fa"), emit: contigs_dup
    tuple val(meta), path("${prefix}/3-mapping/map_to_unitigs.paf.gz")           , emit: mapping
    tuple val(meta), path("${prefix}/*.log")                                     , emit: log
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // Note: MyloAsm works best with FASTQ files for base quality information
    """
    myloasm \\
        $reads \\
        -o ${prefix} \\
        -t $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        myloasm: \$(myloasm --version | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}/alternate_assemblies
    mkdir -p ${prefix}/3-mapping
    touch ${prefix}/assembly_primary.fa
    touch ${prefix}/final_contig_graph.gfa
    touch ${prefix}/alternate_assemblies/assembly_alternate.fa
    touch ${prefix}/alternate_assemblies/duplicated_contigs.fa
    echo | gzip > ${prefix}/3-mapping/map_to_unitigs.paf.gz
    touch ${prefix}/myloasm_1.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        myloasm: \$(myloasm --version | sed 's/.* //')
    END_VERSIONS
    """
}
