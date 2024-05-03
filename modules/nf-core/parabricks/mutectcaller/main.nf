process PARABRICKS_MUTECTCALLER {
    tag "$meta.id"
    label 'process_high'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1"

    /*
    NOTE: Parabricks requires the files to be non-symlinked
    Do not change the stageInMode to soft linked! This is default on Nextflow.
    If you change this setting be careful.
    */
    stageInMode "copy"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index),  path(normal_bam), path(normal_bam_index), path(interval_file)
    tuple val(ref_meta), path(fasta)
    path panel_of_normals
    path panel_of_normals_index

    output:
    tuple val(meta), path("*.vcf.gz"),       emit: vcf
    tuple val(meta), path("*.vcf.gz.stats"), emit: stats
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def prepon_command = panel_of_normals ? "cp -L $panel_of_normals_index `readlink -f $panel_of_normals`.tbi && pbrun prepon --in-pon-file $panel_of_normals" : ""
    def postpon_command = panel_of_normals ? "pbrun postpon --in-vcf ${prefix}.vcf.gz --in-pon-file $panel_of_normals --out-vcf ${prefix}_annotated.vcf.gz" : ""
    """

    # if panel of normals specified, run prepon
    $prepon_command

    pbrun \\
        mutectcaller \\
        --ref $fasta \\
        --in-tumor-bam $tumor_bam \\
        --tumor-name ${meta.tumor_id} \\
        --out-vcf ${prefix}.vcf.gz \\
        $interval_file_command \\
        --num-gpus $task.accelerator.request \\
        $args

    # if panel of normals specified, run postpon
    $postpon_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def postpon_command = panel_of_normals ? "touch ${prefix}_annotated.vcf.gz" : ""
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.stats
    $postpon_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
