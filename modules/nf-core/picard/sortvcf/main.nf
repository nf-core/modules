process PICARD_SORTVCF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b474c6f12c0502377b95062b75ef4b7778864b1353c7fc1d5a3b1f3b3017fd2e/data'
:         'community.wave.seqera.io/library/picard:3.5.0--842d4c70c98af9b4' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(dict)

    output:
    tuple val(meta), path("*_sorted.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('picard'), eval("picard SortVcf --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_dict = dict ? "--SEQUENCE_DICTIONARY ${dict}" : ""
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard SortVcf] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        SortVcf \\
        --INPUT ${vcf} \\
        ${args} \\
        ${seq_dict} \\
        ${reference} \\
        --OUTPUT ${prefix}_sorted.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_sorted.vcf.gz
    touch ${prefix}.bam.bai
    touch ${prefix}.MarkDuplicates.metrics.txt
    """
}
