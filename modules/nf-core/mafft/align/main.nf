process MAFFT_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0':
        'biocontainers/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(add)
    tuple val(meta3), path(addfragments)
    tuple val(meta4), path(addfull)
    tuple val(meta5), path(addprofile)
    tuple val(meta6), path(addlong)
    val(compress)

    output:
    tuple val(meta), path("*.fas{.gz,}"), emit: fas
    tuple val("${task.process}"), val("mafft"), eval("mafft --version 2>&1 | sed 's/ (.*) //g'"), topic: versions, emit: versions_mafft
    tuple val("${task.process}"), val("pigz"), eval("pigz --version 2>&1 | sed 's/pigz //g'")   , topic: versions, emit: versions_pigz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    add          = add             ? "--add <(unpigz -cdf ${add})"                   : ''
    addfragments = addfragments    ? "--addfragments <(unpigz -cdf ${addfragments})" : ''
    addfull      = addfull         ? "--addfull <(unpigz -cdf ${addfull})"           : ''
    addprofile   = addprofile      ? "--addprofile <(unpigz -cdf ${addprofile})"     : ''
    addlong      = addlong         ? "--addlong <(unpigz -cdf ${addlong})"           : ''
    write_output = compress ? " | pigz -cp ${task.cpus} > ${prefix}.fas.gz" : "> ${prefix}.fas"
    // this will not preserve MAFFTs return value, but mafft crashes when it receives a process substitution
    if ("$fasta" == "${prefix}.fas" ) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    mafft \\
        --thread ${task.cpus} \\
        ${add} \\
        ${addfragments} \\
        ${addfull} \\
        ${addprofile} \\
        ${addlong} \\
        ${args} \\
        ${fasta} \\
        ${write_output}
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fas" ) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    echo ${args}

    if [[ "$compress" == "true" ]]; then
        echo "" | pigz -cp ${task.cpus} > ${prefix}.fas.gz
    else
        touch ${prefix}.fas
    fi
    """

}
