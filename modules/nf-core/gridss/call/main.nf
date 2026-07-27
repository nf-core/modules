process GRIDSS_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h50ea8bc_3':
        'quay.io/biocontainers/gridss:2.13.2--h50ea8bc_3' }"

    input:
    tuple val(meta), path(bam), path(bai), path(preprocess_dir), path(assemble_dir), path(assemble_bam)
    tuple val(meta2), path(fasta), path(fasta_fai), path(bwa_index)
    tuple val(meta3), path(gridss_config)

    output:
    tuple val(meta), path("*.sv.gridss.vcf.gz"),  emit: vcf
    tuple val("${task.process}"), val('gridss'), eval("CallVariants --version 2>&1 | sed 's/-gridss\$//'")  , topic: versions, emit: versions_gridss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_config = gridss_config ? "-c ${gridss_config}" : ""
    def bam_list = bam instanceof List ? bam : [bam]

    def index_files = bwa_index instanceof List ? bwa_index : [bwa_index]
	// GRIDSS requires all BWA index files to have the exact same basename as the reference fasta
	def link_cmds = index_files.collect { idx -> "ln -sf ${idx} ./${fasta}.${idx.extension}" }.join('\n')
    """
    ${link_cmds}

    gridss ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --steps call \\
        --reference ${fasta} \\
        --workingdir  "." \\
        --assembly ${assemble_bam} \\
        --output ${prefix}.sv.gridss.vcf.gz \\
        --threads ${task.cpus} ${arg_config} ${bam_list.join(' ')}

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.sv.gridss.vcf.gz
    """
}
