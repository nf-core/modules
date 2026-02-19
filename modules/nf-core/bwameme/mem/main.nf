process BWAMEME_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9ddd41b93c5e182db9d643ca266dd1677e59593a9cb49904b982ff45ad5aa8c3/data':
        'community.wave.seqera.io/library/bwa-meme_mbuffer_samtools:03f3f60b6c289776' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    val   sort_bam
    val   mbuffer
    val   samtools_threads

    output:
    tuple val(meta), path("*.sam")  , emit: sam , optional:true
    tuple val(meta), path("*.bam")  , emit: bam , optional:true
    tuple val(meta), path("*.cram") , emit: cram, optional:true
    tuple val(meta), path("*.crai") , emit: crai, optional:true
    tuple val(meta), path("*.csi")  , emit: csi , optional:true
    tuple val("${task.process}"), val('bwameme'), val('1.0.6'), topic: versions, emit: versions_bwameme
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    if (!mbuffer) {
        log.info '[bwameme-mbuffer] Memory for mbuffer is not set - defaulting to 3GB for mbuffer.'
        mbuffer_mem = 3072
    } else {
        mbuffer_mem = mbuffer
    }
    if (!samtools_threads) {
        log.info 'Number of threads for samtools is not set - defaulting to 2 threads.'
        threads = 2
    } else {
        threads = samtools_threads
    }
    mbuffer_command   = sort_bam ? "| mbuffer -m ${mbuffer_mem}M" : ""
    mem_per_thread    = sort_bam ? "-m "+ (mbuffer_mem/threads).intValue()+"M" : ""
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-meme \\
        mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        $mbuffer_command \\
        | samtools $samtools_command $args2 $mem_per_thread -@ $threads ${reference} -o ${prefix}.${extension} -
    """

    stub:

    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    } else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }
    """
    touch ${prefix}.${extension}
    ${create_index}
    """
}
