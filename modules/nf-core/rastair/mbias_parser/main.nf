process RASTAIR_MBIAS_PARSER {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rastair:0.8.2--r44h4349ce8_2' :
        'quay.io/biocontainers/rastair:0.8.2--r44h4349ce8_2' }"

    input:
    tuple val(meta), path(rastair_mbias_txt)

    output:
    tuple val(meta), path("*.rastair_mbias_processed.pdf"),         emit: mbias_processed_pdf, optional: true
    tuple val(meta), path("*.rastair_mbias_processed.csv"),         emit: mbias_processed_csv
    tuple val(meta), env(trim_OT), env(trim_OB),                    emit: mbias_processed_str
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_mbias.R --pdf -o ${prefix}.rastair_mbias_processed.pdf ${rastair_mbias_txt} > ${prefix}.rastair_mbias_processed.txt

    parse_mbias.R ${prefix}.rastair_mbias_processed.txt ${prefix}.rastair_mbias_processed.csv
    trim_OT=\$(head -n 1 ${prefix}.rastair_mbias_processed.csv)
    trim_OB=\$(head -n 2 ${prefix}.rastair_mbias_processed.csv | tail -n 1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
