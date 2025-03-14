include { input_args } from '../simpleaf/index/main.nf'
process HIFIASM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiasm:0.24.0--h5ca1c30_0' :
        'biocontainers/hifiasm:0.24.0--h5ca1c30_0' }"

    input:
    tuple val(meta) , path(long_reads)        , path(ul_reads)
    tuple val(meta1), path(paternal_kmer_dump), path(maternal_kmer_dump)
    tuple val(meta2), path(hic_read1)         , path(hic_read2)
    tuple val(meta3), path(bin_files)

    output:
    tuple val(meta), path("*.r_utg.gfa")                             , emit: raw_unitigs
    tuple val(meta), path("*.ec.bin")                                , emit: corrected_reads,   optional: true
    tuple val(meta), path("*.ovlp.source.bin")                       , emit: source_overlaps,   optional: true
    tuple val(meta), path("*.ovlp.reverse.bin")                      , emit: reverse_overlaps,  optional: true
    tuple val(meta), path("*.p_utg.gfa")                             , emit: processed_unitigs, optional: true
    tuple val(meta), path("${prefix}.{p_ctg,bp.p_ctg,hic.p_ctg}.gfa"), emit: primary_contigs  , optional: true
    tuple val(meta), path("${prefix}.{a_ctg,hic.a_ctg}.gfa")         , emit: alternate_contigs, optional: true
    tuple val(meta), path("${prefix}.*.hap1.p_ctg.gfa")              , emit: hap1_contigs     , optional: true
    tuple val(meta), path("${prefix}.*.hap2.p_ctg.gfa")              , emit: hap2_contigs     , optional: true
    tuple val(meta), path("${prefix}.stderr.log")                    , emit: log
    path  "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def long_reads_sorted = long_reads instanceof List ? long_reads.sort{ it.name } : long_reads
    def ul_reads_sorted = ul_reads instanceof List ? ul_reads.sort{ it.name } : ul_reads
    def ultralong = ul_reads ? "--ul ${ul_reads_sorted}" : ""

    if([paternal_kmer_dump, maternal_kmer_dump].any() && [hic_read1, hic_read2].any()) {
        log.error("ERROR: hifiasm trio binning mode and Hi-C phasing can not be used at the same time.")
    }

    def input_trio = ""
    if([paternal_kmer_dump, maternal_kmer_dump].any()) {
        if(![paternal_kmer_dump, maternal_kmer_dump].every()) {
            log.error("ERROR: Either the maternal or paternal kmer dump is missing!")
        } else {
            input_trio = "-1 ${paternal_kmer_dump} -2 ${maternal_kmer_dump}"
        }
    }

    def input_hic = ""
    if([hic_read1, hic_read2].any()) {
        if(![hic_read1, hic_read2].every()) {
            log.error("ERROR: Either the forward or reverse Hi-C reads are missing!")
        } else {
            input_hic = "--h1 ${hic_read1} --h2 ${hic_read2}"
        }
    }
    """
    hifiasm \\
        $args \\
        -t ${task.cpus} \\
        ${input_trio} \\
        ${input_hic} \\
        ${ultralong} \\
        -o ${prefix} \\
        ${long_reads_sorted} \\
        2> >( tee ${prefix}.stderr.log >&2 )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: \$(hifiasm --version 2>&1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.r_utg.gfa
    touch ${prefix}.ec.bin
    touch ${prefix}.ovlp.source.bin
    touch ${prefix}.ovlp.reverse.bin
    touch ${prefix}.bp.p_ctg.gfa
    touch ${prefix}.p_utg.gfa
    touch ${prefix}.p_ctg.gfa
    touch ${prefix}.a_ctg.gfa
    touch ${prefix}.hap1.p_ctg.gfa
    touch ${prefix}.hap2.p_ctg.gfa
    touch ${prefix}.stderr.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: \$(hifiasm --version 2>&1)
    END_VERSIONS
    """
}
