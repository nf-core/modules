process PARABRICKS_MUTECTCALLER {
    tag "$meta.id"
    label 'process_high'
    memory "30 GB"
    cpus 16
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
    tuple val(meta), path(fasta), path(tumor_bam)
    tuple val(meta2), path(normal_bam)
    path interval_file

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.stats"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def normal_bam_command = normal_bam ? "--in-normal-bam $normal_bam --normal-name $meta2.id" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    """
    pbrun \\
        mutectcaller \\
        --ref $fasta \\
        --in-tumor-bam $tumor_bam \\
        $normal_bam_command \\
        --tumor-name ${prefix} \\
        --out-vcf ${prefix}.vcf.gz \\
        $interval_file_command \\
        --num-gpus $task.accelerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}

// TODO
// panel of normals features
// additional mutect arguments
