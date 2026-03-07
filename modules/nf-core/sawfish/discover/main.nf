process SAWFISH_DISCOVER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca71c93b472a8b9a7701a744a5f123e4474ce9bf1e0d110aae5b84b5134dd74c/data' :
        'community.wave.seqera.io/library/sawfish:2.2.0--430c21f2b465b4f7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(expected_cn_bed)
    tuple val(meta4), path(maf_vcf)
    tuple val(meta5), path(cnv_exclude_regions_bed)

    output:
    tuple val(meta), path("${prefix}/assembly.regions.bed")         , emit: assembly_regions
    tuple val(meta), path("${prefix}/candidate.sv.bcf")             , emit: candidate_sv_bcf
    tuple val(meta), path("${prefix}/candidate.sv.bcf.csi")         , emit: candidate_sv_bcf_csi
    tuple val(meta), path("${prefix}/contig.alignment.bam")         , emit: contig_alignment_bam
    tuple val(meta), path("${prefix}/contig.alignment.bam.csi")     , emit: contig_alignment_bam_csi
    tuple val(meta), path("${prefix}/copynum.bedgraph")             , emit: copynum_bedgraph         , optional: true
    tuple val(meta), path("${prefix}/copynum.mpack")                , emit: copynum_mpack            , optional: true
    tuple val(meta), path("${prefix}/debug.breakpoint_clusters.bed"), emit: debug_breakpoint_clusters
    tuple val(meta), path("${prefix}/debug.cluster.refinement.txt") , emit: debug_cluster_refinement
    tuple val(meta), path("${prefix}/discover.settings.json")       , emit: discover_settings
    tuple val(meta), path("${prefix}/genome.gclevels.mpack")        , emit: genome_gclevels          , optional: true
    tuple val(meta), path("${prefix}/max.depth.bed")                , emit: max_depth
    tuple val(meta), path("${prefix}/run.stats.json")               , emit: run_stats
    tuple val(meta), path("${prefix}/sample.gcbias.mpack")          , emit: sample_gcbias            , optional: true
    tuple val(meta), path("${prefix}/sawfish.log")                  , emit: log
    tuple val(meta), path("${prefix}/depth.mpack")                  , emit: depth_mpack
    tuple val(meta), path("${prefix}/maf.mpack")                    , emit: maf_mpack                , optional: true
    tuple val(meta), path("${prefix}/expected.copy.number.bed")     , emit: expected_cn              , optional: true
    tuple val(meta), path("${prefix}")                              , emit: discover_dir
    path("versions.yml")                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    def expected_cn = expected_cn_bed ? "--expected-cn ${expected_cn_bed}" : ""
    def maf = maf_vcf ? "--maf ${maf_vcf}" : ""
    def cnv_exclude_regions = cnv_exclude_regions_bed ? "--cnv-excluded-regions ${cnv_exclude_regions_bed}" : ""

    """
    sawfish \\
        discover \\
        --disable-path-canonicalization \\
        --ref $fasta \\
        --bam $bam \\
        --threads $task.cpus \\
        $args \\
        $expected_cn \\
        $cnv_exclude_regions \\
        $maf \\
        --output-dir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sawfish: \$(sawfish --version | sed 's/sawfish //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def expected_cn = expected_cn_bed ? "touch ${prefix}/expected.copy.number.bed": ''
    def maf_mpack = maf_vcf ? "touch  ${prefix}/maf.mpack": ''
    def copynum_files = !args.contains('--disable-cnv') ? """
        touch ${prefix}/copynum.bedgraph
        touch ${prefix}/copynum.mpack
        touch ${prefix}/sample.gcbias.mpack
        touch ${prefix}/genome.gclevels.mpack
    """ : ''

    """
    mkdir -p ${prefix}
    touch ${prefix}/candidate.sv.bcf
    touch ${prefix}/candidate.sv.bcf.csi
    touch ${prefix}/assembly.regions.bed
    touch ${prefix}/contig.alignment.bam
    touch ${prefix}/contig.alignment.bam.csi
    touch ${prefix}/debug.breakpoint_clusters.bed
    touch ${prefix}/debug.cluster.refinement.txt
    touch ${prefix}/discover.settings.json
    touch ${prefix}/max.depth.bed
    touch ${prefix}/run.stats.json
    touch ${prefix}/sawfish.log
    touch ${prefix}/depth.mpack
    touch ${prefix}/copynum.bedgraph
    touch ${prefix}/copynum.mpack
    touch ${prefix}/sample.gcbias.mpack
    touch ${prefix}/genome.gclevels.mpack
    touch ${prefix}/expected.copy.number.bed

    ${copynum_files}
    ${expected_cn}
    ${maf_mpack}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sawfish: \$(sawfish --version | sed 's/sawfish //g')
    END_VERSIONS
    """
}
