process KMA_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fc6c961562aef21c24b4f2330d9cd7e9bbda162b0d584a5cd5428e0b725e0d6/data':
        'community.wave.seqera.io/library/kma:1.5.0--eb093e0381fb59ea' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("kmaindex"), emit: index
    tuple val("${task.process}"), val('kma'), eval('kma_index -v 2>&1 | sed "s/^KMA_index-//"'), emit: versions_kma, topic: versions


    script:
    def prefix  = task.ext.prefix ?: "${fasta.baseName}"
    def args    = task.ext.args ?: ''
    """
    mkdir kmaindex
    kma \\
        index \\
        -i ${fasta} \\
        -o kmaindex/${prefix} \\
        $args

    """

    stub:
    def prefix  = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir kmaindex

    touch kmaindex/${prefix}.comp.b
    touch kmaindex/${prefix}.length.b
    touch kmaindex/${prefix}.name
    touch kmaindex/${prefix}.seq.b
    """
}
