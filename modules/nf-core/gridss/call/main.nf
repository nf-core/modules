process GRIDSS_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h50ea8bc_3':
        'quay.io/biocontainers/gridss:2.13.2--h50ea8bc_3' }"

    input:
    tuple val(meta), path(bams), path(bais), path(preprocess_dir), path(assemble_dir), path(assemble_bam)
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
    def bams_list = bams instanceof List ? bams : [bams]

    """
    # GRIDSS requires all BWA index files to have the exact
    # same basename as the reference fasta. This is a hard
    # requirement of the tool - it will fail to find indices otherwise.
    index_files=(\$(find -L "${bwa_index}" -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)\$'))

    for index_file in "\${index_files[@]}"; do
        ln -sf "\$index_file" "./${fasta}.\${index_file##*.}"
    done

    gridss ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --steps call \\
        --reference ${fasta} \\
        --workingdir  "." \\
        --assembly ${assemble_bam} \\
        --output ${prefix}.sv.gridss.vcf.gz \\
        --threads ${task.cpus} ${arg_config} ${bams_list.join(' ')}

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.sv.gridss.vcf.gz
    """
}
