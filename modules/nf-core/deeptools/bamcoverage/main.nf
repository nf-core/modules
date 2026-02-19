process DEEPTOOLS_BAMCOVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0':
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta) , path(input)   , path(input_index)
    path(fasta)
    path(fasta_fai)
    tuple val(meta2), path(blacklist)

    output:
    tuple val(meta), path("*.bigWig")  , emit: bigwig  , optional: true
    tuple val(meta), path("*.bedgraph"), emit: bedgraph, optional: true
    tuple val("${task.process}"), val('deeptools'), eval('bamCoverage --version | sed "s/bamCoverage //g"') , emit: versions_deeptools, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'") , emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def blacklist_cmd = blacklist ? "--blackListFileName ${blacklist}" : ""
    def extension = args.contains("--outFileFormat bedgraph") || args.contains("-of bedgraph") ? "bedgraph" : "bigWig"

    // cram_input is currently not working with deeptools
    // therefore it's required to convert cram to bam first
    def is_cram = input.Extension == "cram" ? true : false
    def input_out = is_cram ? input.BaseName + ".bam" : "${input}"
    def fai_reference = fasta_fai ? "--fai-reference ${fasta_fai}" : ""

    if (is_cram){
        """
        samtools view -T $fasta $input $fai_reference -@ $task.cpus -o $input_out
        samtools index -b $input_out -@ $task.cpus

        bamCoverage \\
            --bam $input_out \\
            $args \\
            --numberOfProcessors ${task.cpus} \\
            --outFileName ${prefix}.${extension} \\
            $blacklist_cmd
        """
    }
    else {
        """
        bamCoverage \\
            --bam $input_out \\
            $args \\
            --numberOfProcessors ${task.cpus} \\
            --outFileName ${prefix}.${extension} \\
            $blacklist_cmd
        """
    }

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--outFileFormat bedgraph") || args.contains("-of bedgraph") ? "bedgraph" : "bigWig"
    """
    touch ${prefix}.${extension}
    """
}
