process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_high'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"

    input:
    tuple val(meta), path(reads), path(interval_file)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    path known_sites

    output:
    tuple val(meta), path("*.bam")                , emit: bam
    tuple val(meta), path("*.bai")                , emit: bai
    path "versions.yml"                           , emit: versions
    path "qc_metrics", optional:true              , emit: qc_metrics
    path("*.table"), optional:true                , emit: bqsr_table
    path("duplicate-metrics.txt"), optional:true  , emit: duplicate_metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    def known_sites_command = known_sites ? known_sites.collect{"--knownSites $it"}.join(' ') : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    """

    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    # index and fasta need to be in the same dir as regular files (not symlinks)
    #   and have the same base name for pbrun to function
    #   here we copy the index into the staging dir of fasta
    FASTA_PATH=`readlink -f $fasta`
    cp \$INDEX.amb \$FASTA_PATH.amb
    cp \$INDEX.ann \$FASTA_PATH.ann
    cp \$INDEX.bwt \$FASTA_PATH.bwt
    cp \$INDEX.pac \$FASTA_PATH.pac
    cp \$INDEX.sa \$FASTA_PATH.sa

    pbrun \\
        fq2bam \\
        --ref $fasta \\
        $in_fq_command \\
        --read-group-sm $meta.id \\
        --out-bam ${prefix}.bam \\
        $known_sites_command \\
        $known_sites_output \\
        $interval_file_command \\
        --num-gpus $task.accelerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    def known_sites_command = known_sites ? known_sites.collect{"--knownSites $it"}.join(' ') : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def metrics_output_command = args = "--out-duplicate-metrics duplicate-metrics.txt" ? "touch duplicate-metrics.txt" : ""
    def known_sites_output_command = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output_command = args = "--out-qc-metrics-dir qc_metrics " ? "mkdir qc_metrics && touch qc_metrics/alignment.txt" : ""
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    $metrics_output_command
    $known_sites_output_command
    $qc_metrics_output_command
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
