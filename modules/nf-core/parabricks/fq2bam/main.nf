process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    tuple val(meta4), path(interval_file)
    path(known_sites)

    output:
    tuple val(meta), path("*.bam")  , emit: bam              , optional:true
    tuple val(meta), path("*.bai")  , emit: bai              , optional:true
    tuple val(meta), path("*.cram") , emit: cram             , optional:true
    tuple val(meta), path("*.crai") , emit: crai             , optional:true
    tuple val(meta), path("*.table"), emit: bqsr_table       , optional:true
    path("versions.yml")            , emit: versions
    path("*")                       , emit: qc_metrics       , optional:true
    path("*.txt")                   , emit: duplicate_metrics, optional:true

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    def extension = args2.contains("--output-fmt bam") ? "bam" : "cram"
    def known_sites_command = known_sites ? known_sites.collect{"--knownSites $it"}.join(' ') : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    cp $fasta \$INDEX

    pbrun \\
        fq2bam \\
        --ref \$INDEX \\
        $in_fq_command \\
        --out-bam ${prefix}.${extension} \\
        $known_sites_command \\
        $known_sites_output \\
        $interval_file_command \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-fmt bam") ? "bam" : "cram"
    def extension_index = extension == "cram" ? "crai" : "bai"
    def known_sites_output = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    touch ${prefix}.${extension}
    touch ${prefix}.${extension}.${extension_index}
    $known_sites_output
    $qc_metrics_output
    $duplicate_metrics_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
