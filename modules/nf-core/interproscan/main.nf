process INTERPROSCAN {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "quay.io/nf-core/interproscan:5.73-104.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
    //     'biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)
    tuple path(interproscan_db), val(db_version)

    output:
    tuple val(meta), path('*.tsv.gz') , emit: tsv
    tuple val(meta), path('*.xml.gz') , emit: xml
    tuple val(meta), path('*.gff3.gz'), emit: gff3
    tuple val(meta), path('*.json.gz'), emit: json
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -dp (disable precalculation) is on so no online dependency
    """
    interproscan.sh \\
        -cpu $task.cpus \\
        -i $fasta \\
        -dp \\
        ${args} \\
        --output-file-base ${prefix}

    gzip ${prefix}.{tsv,xml,gff3,json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' > ${prefix}.{tsv,xml,gff3,json}
    gzip ${prefix}.{tsv,xml,gff3,json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """
}
