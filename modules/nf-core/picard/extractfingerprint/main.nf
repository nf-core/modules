process PICARD_EXTRACTFINGERPRINT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  haplotype_map
    path  fasta
    path  fasta_fai
    path  sequence_dictionary

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def bam_name = bam.simpleName

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[PICARD ExtractFingerprint] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    picard \\
        -Xmx${avail_mem}M \\
        ExtractFingerprint \\
        --INPUT ${bam} \\
        --HAPLOTYPE_MAP ${haplotype_map} \\
        --OUTPUT ${prefix}.vcf.gz \\
        ${reference} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard --version 2>&1) | sed 's/^.*(PICARD) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    prefix_no_suffix = task.ext.prefix ? prefix.tokenize('.')[0] : "${meta.id}"
    """
    touch ${prefix_no_suffix}.bam
    touch ${prefix_no_suffix}.cram
    touch ${prefix_no_suffix}.cram.crai
    touch ${prefix_no_suffix}.bai
    touch ${prefix}.metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard --version 2>&1) | sed 's/^.*(PICARD) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
