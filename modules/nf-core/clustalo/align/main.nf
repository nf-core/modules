process CLUSTALO_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4cefc38542f86c17596c29b35a059de10387c6a7:adbe4fbad680f9beb083956d79128039a727e7b3-0':
        'quay.io/biocontainers/mulled-v2-4cefc38542f86c17596c29b35a059de10387c6a7:adbe4fbad680f9beb083956d79128039a727e7b3-0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    path  hmm_in
    path  hmm_batch
    path  profile1
    path  profile2
    val   compress

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    tuple val("${task.process}"), val('clustalo'), eval('clustalo --version'), emit: versions_clustalo, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def options_tree = tree ? "--guidetree-in=$tree" : ""
    def fhmm_in      = hmm_in    ? "--hmm-in=${hmm_in}"       : ""
    def fhmm_batch   = hmm_batch ? "--hmm-batch=${hmm_batch}" : ""
    def fprofile1    = profile1  ? "--profile1=${profile1}"   : ""
    def fprofile2    = profile2  ? "--profile2=${profile2}"   : ""
    def write_output = compress ? "--force -o >(pigz -cp ${task.cpus} > ${prefix}.aln.gz) && wait \$!" : "-o ${prefix}.aln"
    // Process substitution keeps alignment data separate from verbose stdout.
    // Its status is not part of clustalo's status: wait for pigz before success.
    // && preserves clustalo failures; wait propagates compression failures.
    // --force permits opening the existing /dev/fd/<id> path.
    """
    clustalo \
        -i ${fasta} \
        $options_tree \
        ${fhmm_in} \
        ${fhmm_batch} \
        ${fprofile1} \
        ${fprofile2} \
        --threads=${task.cpus} \
        $args \
        $write_output
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln${compress ? '.gz' : ''}
    """
}
