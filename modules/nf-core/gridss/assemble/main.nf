process GRIDSS_ASSEMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h50ea8bc_3':
        'quay.io/biocontainers/gridss:2.13.2--h50ea8bc_3' }"

    input:
    tuple val(meta), path(bams), path(bais), path(preprocess_dirs)
    tuple val(meta2), path(fasta), path(fasta_fai), path(bwa_index)
    tuple val(meta3), path(gridss_config)

    output:
    tuple val(meta), path("*.gridss.working"), path("*.sv.assembly.bam"),  emit: assemble_dir
    tuple val("${task.process}"), val('gridss'), eval("CallVariants --version 2>&1 | sed 's/-gridss//'")  , topic: versions, emit: versions_gridss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_config = gridss_config ? "-c ${gridss_config}" : ""

    def bams_list = bams instanceof List ? bams : [bams]
	def index_files = bwa_index instanceof List ? bwa_index : [bwa_index]
	// GRIDSS requires all BWA index files to have the exact same basename as the reference fasta 
	def link_cmds = index_files.collect { idx -> "ln -sf ${idx} ./${fasta}.${idx.extension}" }.join('\n')
	"""
	${link_cmds}

    gridss ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --steps assemble \\
        --reference ${fasta} \\
        --workingdir "." \\
        --assembly ${prefix}.sv.assembly.bam \\
        --threads ${task.cpus} \\
        ${arg_config} \\
        ${bams_list.join(' ')}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}.sv.assembly.bam.gridss.working

    touch ${prefix}.sv.assembly.bam
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.cigar_metrics
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.coverage.blacklist.bed
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.downsampled_0.bed
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.excluded_0.bed
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.idsv_metrics
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.mapq_metrics
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.quality_distribution_metrics
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.subsetCalled_0.bed
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.sv.bam
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.sv.bam.bai
    touch ${prefix}.sv.assembly.bam.gridss.working/${prefix}.sv.assembly.bam.tag_metrics

    """
}
