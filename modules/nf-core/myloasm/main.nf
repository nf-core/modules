process MYLOASM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b752566c7444cbfddc44bfb332078cf36602fcfeb1e57887cee0d5d6195e1923/data'
        : 'community.wave.seqera.io/library/myloasm:0.5.1--1699da7b77a4bbdd'}"

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
    tuple val("${task.process}"), val('myloasm'), eval("myloasm --version | sed 's/.* //'"), emit: versions_myloasm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    myloasm \\
        ${reads} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}/alternate_assemblies
    mkdir -p ${prefix}/3-mapping
    touch ${prefix}/assembly_primary.fa
    touch ${prefix}/final_contig_graph.gfa
    touch ${prefix}/alternate_assemblies/assembly_alternate.fa
    touch ${prefix}/alternate_assemblies/duplicated_contigs.fa
    echo "" | gzip > ${prefix}/3-mapping/map_to_unitigs.paf.gz
    touch ${prefix}/myloasm_1.log
    """
}
