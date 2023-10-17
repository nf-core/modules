process CHECKM2_PREDICT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0':
        'biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta, stageAs: "input_bins/*")
    path(db)

    output:
    tuple val(meta), path("${prefix}")                   , emit: checkm2_output
    tuple val(meta), path("${prefix}/quality_report.tsv"), emit: checkm2_tsv
    path("versions.yml")                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    db = db ?: "CheckM2_database/uniref100.KO.1.dmnd"
    """
    # download database, hard-coded for now (download option broken in current version https://github.com/chklovski/CheckM2/issues/83)
    if [ ! -f ${db} ]; then
        wget https://zenodo.org/records/5571251/files/checkm2_database.tar.gz?download=1 -O checkm2_database.tar.gz
        tar -xzf checkm2_database.tar.gz
        rm checkm2_database.tar.gz
    fi

    checkm2 \\
        predict \\
        --input ${fasta} \\
        --output-directory ${prefix} \\
        --threads ${task.cpus} \\
        --database_path ${db} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/diamond_output ${prefix}/protein_files
    touch ${prefix}/quality_report.tsv ${prefix}/checkm2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
