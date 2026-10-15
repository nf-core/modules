process DNAAPLER_ALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e78ff2b68785061949b030afb1ff73db6a8d21d5d11c032236f1a1ef93477afa/data'
:         'community.wave.seqera.io/library/dnaapler:1.4.0--da1c8ba3b41e40c6' }"

    input:
    tuple val(meta), path(assembly), path(ignore)
    tuple val(meta2), path(custom_db)

    output:
    tuple val(meta), path("${prefix}_reoriented.fasta")                  , emit: fasta
    tuple val(meta), path("${prefix}_reoriented.gfa")                    , emit: gfa    , optional: true
    tuple val(meta), path("${prefix}_all_reorientation_summary.tsv")     , emit: summary
    tuple val(meta), path("${prefix}_MMseqs2_output.txt")                , emit: hits   , optional: true
    tuple val(meta), path("${prefix}.log")                               , emit: log
    tuple val("${task.process}"), val('dnaapler'), eval("dnaapler --version | sed 's/.*version //'"), topic: versions, emit: versions_dnaapler

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    prefix            = task.ext.prefix ?: "${meta.id}"
    // dnaapler only treats --ignore as a file if the value contains '/' or a '.'
    def ignore_arg    = ignore ? "--ignore ./${ignore}" : ''
    def custom_db_arg = custom_db ? "--custom_db ${custom_db}" : ''
    """
    dnaapler all \\
        --input ${assembly} \\
        --output ${prefix}_dnaapler \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${ignore_arg} \\
        ${custom_db_arg} \\
        ${args}

    mv ${prefix}_dnaapler/${prefix}_reoriented.fasta .
    mv ${prefix}_dnaapler/${prefix}_all_reorientation_summary.tsv .
    mv ${prefix}_dnaapler/dnaapler_*.log ${prefix}.log
    # Not produced for a GFA without circular contigs (no MMseqs2 search is run)
    if [ -f ${prefix}_dnaapler/${prefix}_MMseqs2_output.txt ]; then
        mv ${prefix}_dnaapler/${prefix}_MMseqs2_output.txt .
    fi
    if [ -f ${prefix}_dnaapler/${prefix}_reoriented.gfa ]; then
        mv ${prefix}_dnaapler/${prefix}_reoriented.gfa .
    fi
    """

    stub:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def is_gfa   = assembly.name ==~ /.*\.gfa(\.(gz|bz2|xz))?$/
    """
    echo "${args}"

    touch ${prefix}_reoriented.fasta
    ${is_gfa ? "touch ${prefix}_reoriented.gfa" : ''}
    touch ${prefix}_all_reorientation_summary.tsv
    touch ${prefix}_MMseqs2_output.txt
    touch ${prefix}.log
    """
}
