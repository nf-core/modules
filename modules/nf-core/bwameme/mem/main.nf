process BWAMEME_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9ddd41b93c5e182db9d643ca266dd1677e59593a9cb49904b982ff45ad5aa8c3/data':
        'community.wave.seqera.io/library/bwa-meme_mbuffer_samtools:03f3f60b6c289776' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    val   sort_bam

    output:
    tuple val(meta), path("${prefix}.{sam,bam,cram}"), emit: output
    tuple val(meta), path("${prefix}.{csi,crai}")    , emit: index , optional: true
    tuple val("${task.process}"), val('bwameme'), val('1.0.6'), topic: versions, emit: versions_bwameme
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('mbuffer'), eval("mbuffer --version 2>&1 | sed -n 's/mbuffer version //p'") , topic: versions, emit: versions_mbuffer


    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    // ext.args2 controls mbuffer options; inject default -m if not supplied
    def mbuffer_args = args2.contains('-m') ? args2 : "-m 3072M ${args2}".trim()
    def mbuffer_command = sort_bam ? "| mbuffer ${mbuffer_args}" : ""
    // ext.args3 controls samtools options; inject defaults for -@ and -m (sort only) if not supplied
    def samtools_threads_arg = args3.contains('-@') ? '' : '-@ 3'
    def samtools_mem_arg     = (sort_bam && !args3.contains('-m')) ? '-m 1024M' : ''
    def samtools_args        = "${samtools_mem_arg} ${samtools_threads_arg} ${args3}".trim()
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args3 =~ extension_pattern)
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
        | samtools $samtools_command $samtools_args ${reference} -o ${prefix}.${extension} -
    """

    stub:

    def args3 = task.ext.args3 ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args3 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    def create_index = extension == "cram" ? "touch ${prefix}.crai" :
                       extension == "bam"  ? "touch ${prefix}.csi"  : ""
    """
    touch ${prefix}.${extension}
    ${create_index}
    """
}
