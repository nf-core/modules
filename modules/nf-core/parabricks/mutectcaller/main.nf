process PARABRICKS_MUTECTCALLER {
    tag "$meta.tumor_id"
    label 'process_high'
    memory "30 GB"
    cpus 8
    accelerator 1

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used with Parabricks at the moment. Please use a Docker container."
    }
    if (workflow.containerEngine == 'singularity') {
        exit 1, "Singularity containers cannot be used with Parabricks at the moment. Please use a Docker container."
    }

    container "nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1"
    // container "645946264134.dkr.ecr.us-west-2.amazonaws.com/clara-parabricks:4.0.0-1"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(interval_file)
    tuple val(meta2), path(normal_bam), path(normal_bam_index) 
    path fasta
    path panel_of_normals
    path panel_of_normals_index

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.stats"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumor_id}"
    def normal_bam_command = normal_bam ? "--in-normal-bam $normal_bam --normal-name ${meta2.normal_id}" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def copy_tumor_index_command = tumor_bam_index ? "cp -L $tumor_bam_index `readlink -f $tumor_bam`.bai" : ""
    def copy_normal_index_command = normal_bam_index ? "cp -L $normal_bam_index `readlink -f $normal_bam`.bai" : ""
    def pon_command = panel_of_normals ? "--pon $panel_of_normals" : ""
    def pre_pon_command = panel_of_normals ? "cp -L $panel_of_normals_index `readlink -f $panel_of_normals`.tbi && pbrun prepon --in-pon-file $panel_of_normals" : ""
    """
    # parabricks complains when index is not a regular file in the same directory as the bam
    # copy the index to this path. 
    $copy_tumor_index_command
    $copy_normal_index_command

    # generate pon index file, in panel_of_normals is specified
    $pre_pon_command

    pbrun \\
        mutectcaller \\
        --ref $fasta \\
        --in-tumor-bam $tumor_bam \\
        $normal_bam_command \\
        --tumor-name ${prefix} \\
        --out-vcf ${prefix}.vcf.gz \\
        $interval_file_command \\
        $pon_command \\
        --num-gpus $task.accelerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}

// TODO
// * panel of normals features
// * additional mutect arguments
// * some detection or additional help for the fact that the names specified on --tumor-name 
//     and --normal-name MUST be the same as the sample name specified in the readgroups.
