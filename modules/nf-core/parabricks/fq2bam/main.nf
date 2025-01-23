process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    label 'gpu'
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"

    input:
    tuple val(meta), val(read_group), path ( r1_fastq, stageAs: "?/*"), path ( r2_fastq, stageAs: "?/*")
    tuple path(fasta), path(fai)
    tuple val(meta1), path(index)
    path known_sites
    path interval_file

    output:
    tuple val(meta), path("*.bam"), path("*.bai") , emit: bam_bai
    tuple val(meta), path("qc_metrics/*"), optional: true, emit: qc_metrics
    tuple val(meta), path("*.table"), optional: true, emit: bqsr_table
    tuple val(meta), path("*.duplicate-metrics.txt"), optional: true, emit: duplicate_metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"    
    def known_sites_command = known_sites ? (known_sites instanceof List ? known_sites.collect { "--knownSites $it" }.join(' ') : "--knownSites ${known_sites}") : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? (interval_file instanceof List ? interval_file.collect { "--interval-file $it" }.join(' ') : "--interval-file ${interval_file}") : ""
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''

    def readgroups_string = read_group.collect { rg -> "@RG\\tID:${rg.read_group}__${rg.sample}\\tSM:${rg.sample}\\tPL:${rg.platform}\\tLB:${rg.sample}\\tPU:${rg.read_group}" }

    def in_fq_command = meta.single_end 
        ? (r1_fastq instanceof List 
            ? r1_fastq.collect { "--in-se-fq $it" }.join(' ') 
            : "--in-se-fq ${r1_fastq}"
        )
        : (r1_fastq instanceof List && r2_fastq instanceof List && readgroups_string instanceof List
            ? (r1_fastq.indexed().collect { idx, r1 -> "--in-fq $r1 ${r2_fastq[idx]} \"${readgroups_string[idx]}\"" }).join(' ')
            : "--in-fq ${r1_fastq} ${r2_fastq} \"${readgroups_string.join(' ')}\""
        )
 
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    cp $fasta \$INDEX
    
    pbrun \\
        fq2bam \\
        --ref \$INDEX \\
        $in_fq_command \\
        --read-group-sm $meta.id \\
        --out-bam ${prefix}.bam \\
        $num_gpus \\
        $known_sites_command \\
        $known_sites_output \\
        $interval_file_command \\
        --out-qc-metrics-dir qc_metrics \\
        --out-duplicate-metrics ${prefix}.duplicate-metrics.txt \\
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
    def prefix     = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"    
    def known_sites_command = known_sites ? (known_sites instanceof List ? known_sites.collect { "--knownSites $it" }.join(' ') : "--knownSites ${known_sites}") : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? (interval_file instanceof List ? interval_file.collect { "--interval-file $it" }.join(' ') : "--interval-file ${interval_file}") : ""

    def readgroups_string = read_group.collect { rg -> "@RG\\tID:${rg.read_group}__${rg.sample}\\tSM:${rg.sample}\\tPL:${rg.platform}\\tLB:${rg.sample}\\tPU:${rg.read_group}" }

    def in_fq_command = meta.single_end 
        ? (r1_fastq instanceof List 
            ? r1_fastq.collect { "--in-se-fq $it" }.join(' ') 
            : "--in-se-fq ${r1_fastq}"
        )
        : (r1_fastq instanceof List && r2_fastq instanceof List && readgroups_string instanceof List
            ? (r1_fastq.indexed().collect { idx, r1 -> "--in-fq $r1 ${r2_fastq[idx]} \"${readgroups_string[idx]}\"" }).join(' ')
            : "--in-fq ${r1_fastq} ${r2_fastq} \"${readgroups_string.join(' ')}\""
        )

    def metrics_output_command = args = "--out-duplicate-metrics duplicate-metrics.txt" ? "touch duplicate-metrics.txt" : ""
    def known_sites_output_command = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output_command = args = "--out-qc-metrics-dir qc_metrics " ? "mkdir qc_metrics && touch qc_metrics/alignment.txt" : ""
    """

    echo $in_fq_command

    touch run.log
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
