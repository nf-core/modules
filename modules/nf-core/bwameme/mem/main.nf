process BWAMEME_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ed29b84fa94419f5a7bf6a841ddbcb964768825b:139b5e403886ad278b9ad139174967441c1c6ff3-0':
        'biocontainers/mulled-v2-ed29b84fa94419f5a7bf6a841ddbcb964768825b:139b5e403886ad278b9ad139174967441c1c6ff3-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    val   sort_bam

    output:
    tuple val(meta), path("*.sam")  , emit: sam , optional:true
    tuple val(meta), path("*.bam")  , emit: bam , optional:true
    tuple val(meta), path("*.cram") , emit: cram, optional:true
    tuple val(meta), path("*.crai") , emit: crai, optional:true
    tuple val(meta), path("*.csi")  , emit: csi , optional:true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def mbuffer_mem = 3072
    if (!task.memory) {
        log.info '[bwameme-mbuffer] Available memory not known - defaulting to 3GB for mbuffer. Specify process memory requirements to change this.'
    } else {
        mbuffer_mem = (task.memory.mega*0.5).intValue()
    }
    def mbuffer_command   = sort_bam ? "| mbuffer -m ${mbuffer_mem}M" : ""
    def mem_per_thread    = sort_bam ? "-m "+ (mbuffer_mem/task.cpus).intValue()+"M" : ""
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    def VERSION = '1.0.6' // WARN: Version information provided by tool on CLI is incorrect. Please update this string when bumping container versions.
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-meme \\
        mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        $mbuffer_command \\
        | samtools $samtools_command $args2 $mem_per_thread -@ $task.cpus ${reference} -o ${prefix}.${extension} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameme: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    } else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }
    def VERSION = '1.0.6' // WARN: Version information provided by tool on CLI is incorrect. Please update this string when bumping container versions.
    """
    touch ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameme: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
