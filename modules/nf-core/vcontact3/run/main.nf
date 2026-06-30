process VCONTACT3_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularityOptions ?
        'https://depot.galaxyproject.org/singularity/vcontact3:3.1.6--py39h6e9494a_1' :
        'quay.io/biocontainers/vcontact3:3.1.6--py39h6e9494a_1' }"

    input:
    tuple val(meta), path(genomes)

    output:
    tuple val(meta), path("vcontact3_output/"), emit: results
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    vcontact3 run \\
        -i ${genomes.join(' ')} \\
        -o vcontact3_output/ \\
        --threads ${task.cpus} \\
        ${args}

    cat > versions.yml <<-EOF_VERSIONS
        "${task.process}":
            vcontact3: \$( vcontact3 --version 2>&1 | grep -oP 'vcontact3, version \\K[^\\s]+' )
    EOF_VERSIONS
    """

    stub:
    """
    mkdir -p vcontact3_output/
    touch vcontact3_output/clusters.csv
    touch vcontact3_output/merged_df.csv

    cat > versions.yml <<-EOF_VERSIONS
        "${task.process}":
            vcontact3: 3.1.6
    EOF_VERSIONS
    """
}