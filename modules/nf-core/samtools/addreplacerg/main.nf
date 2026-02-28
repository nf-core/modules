process SAMTOOLS_ADDREPLACERG {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(input), path(index), val(read_group)
    tuple val(meta2), path(fasta), path(fai), path(gzi)

    output:
    tuple val(meta), path("${prefix}.bam")      , emit: bam , optional: true
    tuple val(meta), path("${prefix}.cram")     , emit: cram, optional: true
    tuple val(meta), path("${prefix}.sam")      , emit: sam , optional: true
    tuple val(meta), path("${prefix}.bam.bai")  , emit: bai , optional: true
    tuple val(meta), path("${prefix}.bam.csi")  , emit: csi , optional: true
    tuple val(meta), path("${prefix}.cram.crai"), emit: crai, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference = fasta ? "--reference ${fasta}" : ''
    def read_group_arg = read_group ? "-r ${read_group}" : ''
    file_type = args =~ /--output-fmt\s+sam|-O\s+sam/ ? "sam" :
                args =~ /--output-fmt\s+bam|-O\s+bam/ ? "bam" :
                args =~ /--output-fmt\s+cram|-O\s+cram/ ? "cram" :
                input.getExtension()
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        addreplacerg \\
        --threads ${task.cpus-1} \\
        $args \\
        $read_group_arg \\
        $reference \\
        -o ${prefix}.${file_type} \\
        $input
    """

    stub:
    def file_type = input.getExtension()
    prefix        = task.ext.prefix ?: "${meta.id}"
    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    touch ${prefix}.${file_type}
    """
}
