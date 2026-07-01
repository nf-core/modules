process VCONTACT3_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularityOptions ?
        'docker://quay.io/biocontainers/vcontact3:3.1.6--pyhdfd78af_0' :
        'quay.io/biocontainers/vcontact3:3.1.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(genomes)

    output:
    tuple val(meta), path("vcontact3_output/"), emit: results
    tuple val("${task.process}"), val('vcontact3'), val('3.1.6'), emit: versions_vcontact3, topic: versions

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
    """

    stub:
    """
    mkdir -p vcontact3_output/
    touch vcontact3_output/clusters.csv
    touch vcontact3_output/merged_df.csv
    touch vcontact3_output/vContact_contig_id.txt
    """
}
