// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process DECOUPLER {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "conda-forge::decoupler-py=1.6.0"
    container = "ghcr.io/saezlab/publish-packages/decoupler:sha-2f65a0d"


    input:
    tuple val(meta), path(mat)
    path(net)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*__decoupler.tsv"), emit: dc_out
    // TODO nf-core: List additional required output channels/values here
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '{}'
    def methods = task.ext.methods ?: 'None'
    def source = task.ext.source ?: "source"
    def target = task.ext.target ?: "target"
    def weight = task.ext.weight ?: "weight"
    def min_n = task.ext.min_n ?: 5
    def dense = task.ext.verbose ?: 'False'
    def consensus = task.ext.verbose ?: 'False'
    def verbose = task.ext.verbose ?: 'False'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    import decoupler as dc
    import pandas as pd

    mat = pd.read_csv("${mat}", sep="\t", index_col=0)
    net = pd.read_csv("${net}", sep="\t", index_col=0)

    results = dc.decouple(
        mat=mat,
        net=net,
        methods=${methods},
        source="${source}",
        target="${target}",
        weight="${weight}",
        min_n=int(${min_n}),
        dense=${dense},
        consensus=${consensus},
        verbose=${verbose},
        args=${args}
    )

    for result in results:
        results[result].to_csv(result + "__decoupler.tsv", sep="\t")

    ## VERSIONS FILE
    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("\tdecoupler-py: " + dc.__version__ + "\\n")
    """

    // stub: //dryrun, to-do
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    // //               Have a look at the following examples:
    // //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    // //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // """
    // touch ${prefix}.bam

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    // END_VERSIONS
    // """
}
