process PARABRICKS_FQ2BAM {
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

    // container "nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1"
    container "645946264134.dkr.ecr.us-west-2.amazonaws.com/clara-parabricks:4.0.0-1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    path known_sites
    path interval_file
    val markdups
    val qc_metrics

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions
    path "qc_metrics", optional:true, emit: qc_metrics
    path("*.table"), optional:true, emit: bqsr_table
    path("*-duplicate-metrics.txt"), optional:true, emit: duplicate_metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    def mark_duplicates_command = markdups ? "" : "--no-markdups"
    def mark_duplicates_output = markdups ? "--out-duplicate-metrics ${prefix}-duplicate-metrics.txt" : ""
    def qc_metrics = qc_metrics ? "--out-qc-metrics-dir qc_metrics" : ""
    def known_sites_command = known_sites ? known_sites.collect{"--knownSites $it"}.join(' ') : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    """

    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    # index and fasta need to be in the same dir as files and not symlinks
    # and have the same base name for pbrun to function
    cp \$INDEX.amb ${fasta}.amb 
    cp \$INDEX.ann ${fasta}.ann 
    cp \$INDEX.bwt ${fasta}.bwt 
    cp \$INDEX.pac ${fasta}.pac 
    cp \$INDEX.sa ${fasta}.sa 
    cp -L $fasta ${fasta}2
    mv ${fasta}2 $fasta

    pbrun \\
        fq2bam \\
        --ref $fasta \\
        $in_fq_command \\
        --read-group-sm $meta.id \\
        --out-bam ${prefix}.bam \\
        $mark_duplicates_command \\
        $mark_duplicates_output \\
        $qc_metrics \\
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
}
