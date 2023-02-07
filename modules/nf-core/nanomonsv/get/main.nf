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
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process NANOMONSV_GET {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::nanomonsv=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanomonsv:0.5.0--pyhdfd78af_0':
        'quay.io/biocontainers/nanomonsv:0.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(parse_intermediates), path(control_bam), path(control_bai), path(control_parse_intermediates)
    path reference_fasta

    output:
    tuple val(meta), path("${prefix}.nanomonsv.result.vcf")         , emit: vcf
    tuple val(meta), path("${prefix}.nanomonsv.result.txt")         , emit: text_report
    tuple val(meta), path("${prefix}.nanomonsv.supporting_read.txt"), emit: supp_read
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def control_prefix = task.ext.control_prefix ?: "${meta.id}"
    if (control_prefix) {
        args += " --control_prefix ${control_prefix} --control_bam ${control_bam}"
    }
    """
    nanomonsv get ${args} ${prefix} \
        ${bam} ${reference_fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanomonsv: \$(echo \$(nanomonsv --version 2>&1) | sed 's/^nanomonsv //')
        mafft: \$(echo \$(mafft --version 2>&1) | sed 's/^v//; s/ (.*//')
        racon: \$(echo \$(racon --version 2>&1) | sed 's/^v//')
        tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^tabix (htslib) //; s/ Copyright.*//')
        bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^bgzip (htslib) //; s/ Copyright.*//')
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
