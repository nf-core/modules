// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

process PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pb-cpg-tools:3.0.0--h9ee0642_0':
        'biocontainers/pb-cpg-tools:3.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bed.gz"),      emit: bed
    tuple val(meta), path("*.bed.gz.tbi"),  emit: bed_index
    tuple val(meta), path("*.bw"),          emit: bigwig
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    aligned_bam_to_cpg_scores \\
        --bam ${bam} \\
        --output-prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbcpgtools: \$(aligned_bam_to_cpg_scores --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    echo "" | gzip > ${prefix}.combined.bed.gz
    touch ${prefix}.combined.bed.gz.tbi
    touch ${prefix}.combined.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbcpgtools: \$(aligned_bam_to_cpg_scores --version)
    END_VERSIONS
    """
}
