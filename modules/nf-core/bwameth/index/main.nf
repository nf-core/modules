process BWAMETH_INDEX {
    tag "${fasta}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b8/b806996da500b3cc610c48536ce3a5151a0ff512228e313439534ba08a6acda2/data'
:         'community.wave.seqera.io/library/bwameth:0.2.9--63cb7f0b9f8eae9d' }"

    input:
    tuple val(meta), path(fasta, name: "BwamethIndex/")
    val use_mem2

    output:
    tuple val(meta), path("BwamethIndex"), emit: index
    tuple val("${task.process}"), val('bwameth'), eval("bwameth.py --version | cut -f2 -d' '"), emit: versions_bwameth_index, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_cmd = use_mem2 ? "index-mem2" : "index"
    """
    bwameth.py ${index_cmd} ${fasta}

    rm $fasta
    """

    stub:
    """
    rm ${fasta}

    mkdir -p BwamethIndex/
    touch BwamethIndex/genome.fasta.bwameth.c2t
    touch BwamethIndex/genome.fasta.bwameth.c2t.amb
    touch BwamethIndex/genome.fasta.bwameth.c2t.ann
    touch BwamethIndex/genome.fasta.bwameth.c2t.bwt
    touch BwamethIndex/genome.fasta.bwameth.c2t.pac
    touch BwamethIndex/genome.fasta.bwameth.c2t.sa
    """
}
