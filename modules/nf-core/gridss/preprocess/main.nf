process GRIDSS_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h50ea8bc_3':
        'quay.io/biocontainers/gridss:2.13.2--h50ea8bc_3' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fasta_fai), path(bwa_index)

    output:
    tuple val(meta), path("*.gridss.working"), emit: preprocess_dir
    tuple val("${task.process}"), val('gridss'), eval("CallVariants --version 2>&1 | sed 's/-gridss//'"), topic: versions, emit: versions_gridss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ln -s \$(find -L ${bwa_index} -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)') ./

    gridss \\
        --threads ${task.cpus} \\
        --steps preprocess \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        --reference ${fasta} \\
        ${args} \\
        ${bam}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}.gridss.working/

    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.cigar_metrics
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.computesamtags.changes.tsv
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.coverage.blacklist.bed
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.idsv_metrics
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.insert_size_histogram.pdf
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.insert_size_metrics
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.mapq_metrics
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.sv.bam
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.sv.bam.csi
    touch ${prefix}.gridss.working/${prefix}.gridss.targeted.bam.tag_metrics
    """
}
