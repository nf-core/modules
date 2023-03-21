process PICARD_SORTVCF {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)
    path reference
    path sequence_dict

    output:
    tuple val(meta), path("*_sorted.vcf.gz"), emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_dict = sequence_dict ? "--SEQUENCE_DICTIONARY $sequence_dict" : ""
    def reference = reference ? "--REFERENCE_SEQUENCE $reference" : ""
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
    touch ${prefix}_sorted.vcf.gz
    touch ${prefix}.bam.bai
    touch ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SortVcf --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
