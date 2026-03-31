process DISAMBIGUATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngs-disambiguate:2018.05.03--h06902ac_10':
        'biocontainers/ngs-disambiguate:2018.05.03--h06902ac_10' }"

    input:
    tuple val(meta), path(bam_a), path(bam_b)

    output:
    tuple val(meta), path("*disambiguatedSpeciesA.bam"),                                                      emit: disambiguated_bam_a
    tuple val(meta), path("*disambiguatedSpeciesB.bam"),                                                      emit: disambiguated_bam_b
    tuple val(meta), path("*ambiguousSpeciesA.bam"),                                                          emit: ambiguous_bam_a
    tuple val(meta), path("*ambiguousSpeciesB.bam"),                                                          emit: ambiguous_bam_b
    tuple val(meta), path("*_summary.txt"),                                                                   emit: summary
    tuple val("${task.process}"), val('ngs-disambiguate'), val('2018.05.03'), topic: versions, emit: versions_ngs_disambiguate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ngs_disambiguate \\
        ${bam_a} \\
        ${bam_b} \\
        $args \\
        -s ${prefix} \\
        -o .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.disambiguatedSpeciesA.bam
    touch ${prefix}.disambiguatedSpeciesB.bam
    touch ${prefix}.ambiguousSpeciesA.bam
    touch ${prefix}.ambiguousSpeciesB.bam
    touch ${prefix}_summary.txt
    """
}
