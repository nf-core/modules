process PICARD_SORTVCF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0' :
        'biocontainers/picard:3.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(dict)

    output:
    tuple val(meta), path("*_sorted.vcf.gz"), emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_dict = dict ? "--SEQUENCE_DICTIONARY $dict" : ""
    def reference = fasta ? "--REFERENCE_SEQUENCE $fasta" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard SortVcf] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    picard \\
        SortVcf \\
        -Xmx${avail_mem}M \\
        --INPUT $vcf \\
        $args \\
        $seq_dict \\
        $reference \\
        --OUTPUT ${prefix}_sorted.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SortVcf --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_sorted.vcf.gz
    touch ${prefix}.bam.bai
    touch ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SortVcf --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
