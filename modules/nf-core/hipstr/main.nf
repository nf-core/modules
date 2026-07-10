process HIPSTR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hipstr:0.7--hcf09f9e_0':
        'quay.io/biocontainers/hipstr:0.7--hcf09f9e_0' }"

    input:
    tuple val(meta),  path(alignment_files), path(alignment_indices), path(alignment_file_list)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(ref_regions)
    tuple val(meta4), path(ref_vcf), path(ref_vcf_tbi)
    tuple val(meta5), path(snp_vcf), path(snp_vcf_tbi)
    path stutter_in
    path fam
    path hap_chr_file
    val save_log
    val save_viz
    val save_stutter

    output:
    tuple val(meta), path("*.vcf.gz"),             emit: vcf
    tuple val(meta), path("*.log.txt"),            emit: log, optional: true
    tuple val(meta), path("*.viz.gz"),             emit: viz, optional: true
    tuple val(meta), path("*.stutter_models.txt"), emit: stutter, optional: true
    tuple val("${task.process}"), val('hipstr'), eval("HipSTR --version 2>&1 | sed 's/^HipSTR version //'"), topic: versions, emit: versions_hipstr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (alignment_file_list) {
        log.warn "alignment_file_list is less portable: paths inside it are not staged by Nextflow. Use alignment_files + alignment_indices instead."
    }

    def log_arg = save_log ? "--log ${prefix}.log.txt" : ''
    def viz_arg = save_viz ? "--viz-out ${prefix}.viz.gz" : ''
    def stutter_out_arg = save_stutter ? "--stutter-out ${prefix}.stutter_models.txt" : ''

    def ref_vcf_arg = ref_vcf ? "--ref-vcf $ref_vcf" : ''
    def snp_vcf_arg = snp_vcf ? "--snp-vcf $snp_vcf" : ''
    def stutter_in_arg = stutter_in ? "--stutter-in $stutter_in" : ''
    def fam_arg = fam ? "--fam $fam" : ''
    def hap_chr_file_arg = hap_chr_file ? "--hap-chr-file $hap_chr_file" : ''

    def alignment_file_list_arg = alignment_file_list ? "--bam-files $alignment_file_list" : ''
    def bams_arg = alignment_files ? "--bams ${alignment_files.join(',')}" : ''

    """
    HipSTR \\
        $bams_arg \\
        $alignment_file_list_arg \\
        --fasta $fasta \\
        --regions $ref_regions \\
        --str-vcf ${prefix}.vcf.gz \\
        $ref_vcf_arg \\
        $snp_vcf_arg \\
        $stutter_in_arg \\
        $fam_arg \\
        $hap_chr_file_arg \\
        $log_arg \\
        $viz_arg \\
        $stutter_out_arg \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.log.txt
    echo "" | gzip >  ${prefix}.viz.gz
    touch ${prefix}.stutter_models.txt
    """
}
