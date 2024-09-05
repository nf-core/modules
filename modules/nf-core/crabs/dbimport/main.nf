process CRABS_DBIMPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cutadapt_muscle_vsearch_wget_pruned:0cd5cb1e549e5033':
        'community.wave.seqera.io/library/cutadapt_muscle_vsearch_wget_pruned:04f6c0370c0226c5' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name    = fasta.name.replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    crabs db_import \\
        --input ${fasta_name} \\
        --output ${prefix}.crabsdb.fa \\
        $args

    rm ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs --version | sed -e 's/crabs v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs --version | sed -e 's/crabs v//g')
    END_VERSIONS
    """
}
