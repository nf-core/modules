process GATK4_VARIANTRECALIBRATOR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/3498e3d67f977e18255b0b71560dd548885309f62b736dd2c782b131cf931008/data'
        : 'community.wave.seqera.io/library/gatk4-lite_r-base_r-ggplot2_r-gplots_r-gsalib:aad8b6b5e895ebef'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path resource_vcf
    path resource_tbi
    val labels
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.recal"), emit: recal
    tuple val(meta), path("*.idx"), emit: idx
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R"), emit: plots, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_command = fasta ? "--reference ${fasta} " : ''
    def labels_command = labels.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK VariantRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantRecalibrator \\
        --variant ${vcf} \\
        --output ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        ${reference_command} \\
        --tmp-dir . \\
        ${labels_command} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal
    touch ${prefix}.idx
    touch ${prefix}.tranches
    touch ${prefix}plots.R
    """
}
