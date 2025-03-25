process TRNASCANSE {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trnascan-se:2.0.12--pl5321h7b50bb2_2':
        'biocontainers/trnascan-se:2.0.12--pl5321h7b50bb2_2' }"

    input:
    tuple val(meta), path(fasta)
    val(write_fasta)
    val(write_gff)
    val(write_bed)

    output:
    tuple val(meta), path("*.tsv")   , emit: tsv
    tuple val(meta), path("*.log")   , emit: log
    tuple val(meta), path("*.stats") , emit: stats
    tuple val(meta), path("*.fasta") , emit: fasta , optional: true
    tuple val(meta), path("*.gff")   , emit: gff   , optional: true
    tuple val(meta), path("*.bed")   , emit: bed   , optional: true
    path("versions.yml")             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def fasta_out = write_fasta     ? "-a ${prefix}.fasta" : ""
    def gff       = write_gff       ? "-j ${prefix}.gff" : ""
    def bed       = write_bed       ? "-a ${prefix}.bed" : ""
    def input     = fasta.toString() - ~/\.gz$/
    def unzip     = fasta.getExtension() == "gz" ? "gunzip -c ${fasta} > ${input}" : ""
    def cleanup   = fasta.getExtension() == "gz" ? "rm ${input}" : ""
    """
    ${unzip}

    ## larger genomes can fill up the limited temp space in the singularity container
    ## expected location of the default config file is with the exectuable
    ## copy this and modify to use the working dir as the temp directory
    conf=\$(which tRNAscan-SE).conf
    cp \${conf} trnascan.conf
    sed -i s#/tmp#.#g trnascan.conf

    tRNAscan-SE \\
        --thread ${task.cpus} \\
        -c trnascan.conf \\
        ${args} \\
        -o ${prefix}.tsv \\
        -l ${prefix}.log \\
        -m ${prefix}.stats \\
        ${fasta_out} \\
        ${gff} \\
        ${bed} \\
        ${input}

    if [ -f ${prefix}.fasta ]; then
        ${prefix}.fasta
    fi

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: \$(tRNAscan-SE 2>&1 >/dev/null | awk 'NR==2 {print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def touch_fasta = write_fasta ? "echo '' | gzip > ${prefix}.fasta.gz" : ""
    def touch_gff   = write_gff ? "touch ${prefix}.gff" : ""
    def touch_bed   = write_bed ? "touch ${prefix}.bed" : ""
    """
    touch ${prefix}.tsv
    touch ${prefix}.log
    touch ${prefix}.stats
    ${touch_fasta}
    ${touch_gff}
    ${touch_bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: \$(tRNAscan-SE 2>&1 >/dev/null | awk 'NR==2 {print \$2}')
    END_VERSIONS
    """
}
