process VSEARCH_DEREPLICATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.28.1--h6a68c12_1':
        'biocontainers/vsearch:2.28.1--h6a68c12_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.derep.fasta')   , emit: fasta
    tuple val(meta), path('*.derep.uc')      , emit: clustering
    path "*.derep.log"                       , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vsearch \\
        --derep_fulllength ${fasta} \\
        $args \\
        --relabel "${prefix}." \\
        --uc ${prefix}.derep.uc \\
        --output ${prefix}.derep.fasta 2>&1 | tee ${prefix}.derep.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.derep.fasta
    touch ${prefix}.derep.uc
    touch myfile.derep.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
