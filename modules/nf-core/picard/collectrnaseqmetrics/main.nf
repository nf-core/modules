// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process PICARD_COLLECTRNASEQMETRICS {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::picard=2.27.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    tuple val(meta), path(bam)
    path ref_flat
    path rrna_intervals
    path fasta

    output:
    tuple val(meta), path("*.rna_metrics")  , emit: metrics 
    tuple val(meta), path("*.pdf")          , emit: pdf, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def rrna = rrna_intervals ? "--RIBOSOMAL_INTERVALS ${rrna_intervals}" : ""
    def strand = ""
    if ( meta.strandedness == "forward" || meta.single_end ) {
      strand = "--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND"
    } else if (meta.strandedness == "forward") {
      strand = "--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND"
    } else { strand = "--STRAND_SPECIFICITY NONE" }
    """
    picard \\
        -Xmx${avail_mem}g \\
        CollectRnaSeqMetrics \\
        $args \\
        $reference \\
        $rrna \\
        $strand \\
        --REF_FLAT $ref_flat \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.rna_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectRnaSeqMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rna_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectRnaSeqMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
