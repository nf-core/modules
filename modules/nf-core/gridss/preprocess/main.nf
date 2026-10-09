// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GRIDSS_PREPROCESS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h50ea8bc_3'
        : 'quay.io/biocontainers/gridss:2.13.2--h50ea8bc_3'}"

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

    """
    # GRIDSS requires all BWA index files to have the exact
    # same basename as the reference fasta. This is a hard
    # requirement of the tool - it will fail to find indices otherwise.
    index_files=(\$(find -L "${bwa_index}" -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)\$'))

    for index_file in "\${index_files[@]}"; do
        ln -sf "\$index_file" "./${fasta}.\${index_file##*.}"
    done

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
    def bam_name = bam.getBaseName()

    """
    mkdir -p ${bam_name}.gridss.working/
    cd ${bam_name}.gridss.working/

    touch ${bam_name}.gridss.targeted.bam.cigar_metrics
    touch ${bam_name}.gridss.targeted.bam.computesamtags.changes.tsv
    touch ${bam_name}.gridss.targeted.bam.coverage.blacklist.bed
    touch ${bam_name}.gridss.targeted.bam.idsv_metrics
    touch ${bam_name}.gridss.targeted.bam.insert_size_histogram.pdf
    touch ${bam_name}.gridss.targeted.bam.insert_size_metrics
    touch ${bam_name}.gridss.targeted.bam.mapq_metrics
    touch ${bam_name}.gridss.targeted.bam.sv.bam
    touch ${bam_name}.gridss.targeted.bam.sv.bam.csi
    touch ${bam_name}.gridss.targeted.bam.tag_metrics
    """
}
