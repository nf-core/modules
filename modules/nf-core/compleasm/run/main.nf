process COMPLEASM_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/9392dd488f3f8f78d5567bccb09a0dcf43e9b9ce301321014ba236f7a3c533a2/data'
        : 'community.wave.seqera.io/library/compleasm:0.2.9--52aa1279fc8cbd91'}"

    input:
    tuple val(meta), path(fasta)
    // Required:    lineage to check against, e.g. "bacteria", or "auto" for enabling auto-lineage
    val lineage
    // Recommended: directory holding downloaded BUSCO lineages - downloads if not set
    path library_path

    output:
    tuple val(meta), path("*/summary.txt"), emit: summary
    tuple val(meta), path("*/*/full_table.tsv"), emit: full_table
    tuple val(meta), path("*/*/full_table_busco_format.tsv"), emit: full_table_busco_format
    tuple val(meta), path("*/*/detected_genes.gff"), emit: gff
    tuple val(meta), path("*/*/gene_marker.fasta"), emit: gene_marker
    tuple val(meta), path("*/*/translated_protein.fasta"), emit: translated_protein
    tuple val(meta), path("mb_downloads"), emit: library, optional: true
    tuple val("${task.process}"), val('compleasm'), eval("compleasm --version | sed 's/compleasm //'"), emit: versions_compleasm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def compleasm_lineage = lineage == 'auto' ? '--autolineage' : "--lineage ${lineage}"
    def compleasm_library = library_path ? "-L ${library_path}" : ''
    """
    compleasm run \\
        --threads ${task.cpus} \\
        --assembly_path ${fasta} \\
        --output_dir ${prefix} \\
        ${compleasm_lineage} \\
        ${compleasm_library} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    """
    mkdir -p ${prefix}/${lineage}
    touch ${prefix}/summary.txt
    touch ${prefix}/${lineage}/full_table.tsv
    touch ${prefix}/${lineage}/full_table_busco_format.tsv
    touch ${prefix}/${lineage}/detected_genes.gff
    touch ${prefix}/${lineage}/gene_marker.fasta
    touch ${prefix}/${lineage}/translated_protein.fasta
    """
}
