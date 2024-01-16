process MAFFT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6':
        'biocontainers/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6' }"

    input:
    tuple val(meta),  path(fasta)
    tuple val(meta2), path(add)
    tuple val(meta3), path(addfragments)
    tuple val(meta4), path(addfull)
    tuple val(meta5), path(addprofile)
    tuple val(meta6), path(addlong)

    output:
    tuple val(meta), path("*.fas.gz"), emit: fas
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def add          = add             ? "--add ${add}"                   : ''
    def addfragments = addfragments    ? "--addfragments ${addfragments}" : ''
    def addfull      = addfull         ? "--addfull ${addfull}"           : ''
    def addprofile   = addprofile      ? "--addprofile ${addprofile}"     : ''
    def addlong      = addlong         ? "--addlong ${addlong}"           : ''
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
        | pigz -p ${task.cpus} -c > ${prefix}.fas.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def add          = add             ? "--add ${add}"                   : ''
    def addfragments = addfragments    ? "--addfragments ${addfragments}" : ''
    def addfull      = addfull         ? "--addfull ${addfull}"           : ''
    def addprofile   = addprofile      ? "--addprofile ${addprofile}"     : ''
    def addlong      = addlong         ? "--addlong ${addlong}"           : ''
    if ("$fasta" == "${prefix}.fas" )  error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch ${prefix}.fas.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

}
