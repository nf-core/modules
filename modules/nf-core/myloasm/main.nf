process MYLOASM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/22/225ce5652f6d4f9cb0bc9acb813f87d32abe42fcf0b9a41125c8bb6a84ca3709/data'
        : 'community.wave.seqera.io/library/myloasm:0.7.0--497d9a84b31cc548'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val(meta), path("${prefix}/assembly_primary.fa.gz"), emit: contigs
    tuple val(meta), path("${prefix}/final_contig_graph.gfa.gz"), emit: gfa
    tuple val(meta), path("${prefix}/alternate_assemblies/assembly_alternate.fa.gz"), emit: contigs_alt
    tuple val(meta), path("${prefix}/alternate_assemblies/duplicated_contigs.fa.gz"), emit: contigs_dup
    tuple val(meta), path("${prefix}/3-mapping/map_to_unitigs.paf.gz"), emit: mapping
    tuple val(meta), path("${prefix}/*.log"), emit: log
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

    find ${prefix}/ -name "*.fa" -exec gzip {} \\;
    find ${prefix}/ -name "*.gfa" -exec gzip {} \\;
    find ${prefix}/ -name "*.edges" -exec gzip {} \\;
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}/alternate_assemblies
    mkdir -p ${prefix}/3-mapping
    echo "" | gzip > ${prefix}/assembly_primary.fa.gz
    echo "" | gzip > ${prefix}/final_contig_graph.gfa.gz
    echo "" | gzip > ${prefix}/alternate_assemblies/assembly_alternate.fa.gz
    echo "" | gzip > ${prefix}/alternate_assemblies/duplicated_contigs.fa.gz
    echo "" | gzip > ${prefix}/3-mapping/map_to_unitigs.paf.gz
    touch ${prefix}/myloasm_1.log
    """
}
