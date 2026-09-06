process PANAROO_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/panaroo:1.6.0--pyhdfd78af_0':
        'quay.io/biocontainers/panaroo:1.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${prefix}/*")                      , emit: results
    tuple val(meta), path("${prefix}/core_gene_alignment.aln"), emit: aln, optional: true
    tuple val("${task.process}"), val('panaroo'), eval("panaroo --version 2>&1 | sed 's/^.*panaroo //'") , emit: versions_panaroo, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    panaroo \\
        ${args} \\
        -t ${task.cpus} \\
        -o ${prefix} \\
        -i ${gff}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def core_cmd = args.contains("-a core")

    """
    mkdir -p ${prefix}
    touch ${prefix}/combined_{DNA,protein}_CDS.fasta
    touch ${prefix}/combined_protein_cdhit_out.txt{,.clstr}
    touch ${prefix}/final_graph.gml
    touch ${prefix}/gene_data.csv
    touch ${prefix}/gene_presence_absence.{Rtab,csv}
    touch ${prefix}/gene_presence_absence_roary.csv
    touch ${prefix}/pan_genome_reference.fa
    touch ${prefix}/pre_filt_graph.gml
    touch ${prefix}/struct_presence_absence.Rtab
    touch ${prefix}/summary_statistics.txt

    if [ "${core_cmd}" == "true" ]; then
        touch ${prefix}/aligned_gene_sequences
        touch ${prefix}/alignment_entropy.csv
        touch ${prefix}/core_alignment{,_filtered}_header.embl
        touch ${prefix}/core_gene_alignment{,_filtered}.aln
    fi
    """
}
