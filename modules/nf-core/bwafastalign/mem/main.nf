process BWAFASTALIGN_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f8c975324a12014c8a817c2c1ad0cd68b077cf09c4370717589177b262dcd1dc/data':
        'community.wave.seqera.io/library/bwa-fastalign_mbuffer_samtools:35f24ce8addcd26b'}"

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
    tuple val("${task.process}"), val('bwafastalign'), val('1.0.0'), topic: versions, emit: versions_bwafastalign
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed '1!d;s/.* //'") , topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('mbuffer'), eval("mbuffer --version 2>&1 | sed -n 's/mbuffer //p'") , topic: versions, emit: versions_mbuffer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    if (!mbuffer) {
        log.info '[bwafastalign-mbuffer] Memory for mbuffer is not set - defaulting to 3GB for mbuffer.'
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

    bwa-fastalign \\
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
