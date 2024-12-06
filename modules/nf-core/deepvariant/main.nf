def deprecation_message = """
WARNING: The deepvariant process has been moved into deepvariant/rundeepvariant.

Reason:
A subworkflow "deepvariant" was added, to split DeepVariant into three processes, to help optimise
the usage of GPU resources (https://github.com/nf-core/modules/pull/6172). The subworkflow can be
used instead of this module. Alternatively, it is possible to use the subcommand "rundeepvariant" in
this module. "rundeepvariant" (the process DEEPVARIANT_RUNDEEPVARIANT) is the exact same as this
top-level process (DEEPVARIANT) used to be.

The processing stages used by the subworkflow are implemented as module subcommands, e.g. makeexamples.
"""


process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/deepvariant:1.8.0"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)
    tuple val(meta5), path(par_bed)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    ,  emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"),  emit: gvcf_tbi
    path "versions.yml"                            ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assert false: deprecation_message

}
