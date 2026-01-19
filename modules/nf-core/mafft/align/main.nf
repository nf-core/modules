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
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def add_opt          = add             ? "--add <(unpigz -cdf ${add})"                   : ''
    def addfragments_opt = addfragments    ? "--addfragments <(unpigz -cdf ${addfragments})" : ''
    def addfull_opt      = addfull         ? "--addfull <(unpigz -cdf ${addfull})"           : ''
    def addprofile_opt   = addprofile      ? "--addprofile <(unpigz -cdf ${addprofile})"     : ''
    def addlong_opt      = addlong         ? "--addlong <(unpigz -cdf ${addlong})"           : ''
    def write_output = compress ? " | pigz -cp ${task.cpus} > ${prefix}.fas.gz" : "> ${prefix}.fas"
    // this will not preserve MAFFTs return value, but mafft crashes when it receives a process substitution
    if ("$fasta" == "${prefix}.fas" ) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    mafft \\
        --thread ${task.cpus} \\
        ${add_opt} \\
        ${addfragments_opt} \\
        ${addfull_opt} \\
        ${addprofile_opt} \\
        ${addlong_opt} \\
        ${args} \\
        ${fasta} \\
        ${write_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fas" )  error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch ${prefix}.fas${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

}
