// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process ARCHR_CREATEARROWFILES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "archr" // TODO update

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path("*fragments.tsv.gz")

    output:
    tuple val(meta), path("*.arrow")                , emit: arrow file
    tuple val(meta), path("*.rds")                  , emit: arrow metadata
    tuple val(meta), path(".pdf")                   , emit: quality control plots
    tuple val(meta), path("*.R_sessionInfo.log")    , emit: session_info
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    """
    template createarrowfiles.r
    """

    stub:
    """
    touch QualityControl/${meta.id}/${meta.id}-Pre-Filter-Metadata.rds
    touch QualityControl/${meta.id}/${meta.id}-Fragment_Size_Distribution.pdf
    touch QualityControl/${meta.id}/${meta.id}-TSS_by_Unique_Frags.pdf
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ArchR: \$(Rscript -e "library(ArchR); cat(as.character(packageVersion('ArchR')))")
    END_VERSIONS
    """
}
