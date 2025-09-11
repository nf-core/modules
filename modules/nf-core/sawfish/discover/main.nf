process SAWFISH_DISCOVER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/72/72b081ac73287f98b0453a7e2fbd66e581f9eed68fe91289b9e4189e639fa6d9/data' :
        'community.wave.seqera.io/library/sawfish:2.0.5--422cc4cf3cd63e02' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(expected_cn_bed)
    tuple val(meta4), path(maf_vcf)
    tuple val(meta5), path(cnv_exclude_regions)

    output:
    tuple val(meta), path("versions.yml")                                  , emit: versions
    tuple val(meta), path("*/assembly.regions.bed")           , emit: assembly_regions
    tuple val(meta), path("*/candidate.sv.bcf")               , emit: candidate_sv_bcf
    tuple val(meta), path("*/candidate.sv.bcf.csi")           , emit: candidate_sv_bcf_csi
    tuple val(meta), path("*/contig.alignment.bam")           , emit: contig_alignment_bam
    tuple val(meta), path("*/contig.alignment.bam.csi")       , emit: contig_alignment_bam_csi
    tuple val(meta), path("*/copynum.bedgraph")               , emit: copynum_bedgraph          , optional: true
    tuple val(meta), path("*/copynum.mpack")                  , emit: copynum_mpack             , optional: true
    tuple val(meta), path("*/debug.breakpoint_clusters.bed")  , emit: debug_breakpoint_clusters
    tuple val(meta), path("*/debug.cluster.refinement.txt")   , emit: debug_cluster_refinement
    tuple val(meta), path("*/discover.settings.json")         , emit: discover_settings
    tuple val(meta), path("*/genome.gclevels.mpack")          , emit: genome_gclevels
    tuple val(meta), path("*/max.depth.bed")                  , emit: max_depth
    tuple val(meta), path("*/run.stats.json")                 , emit: run_stats
    tuple val(meta), path("*/sample.gcbias.mpack")            , emit: sample_gcbias             , optional: true
    tuple val(meta), path("*/sawfish.log")                    , emit: log
    tuple val(meta), path("*/depth.mpack")                    , emit: depth_mpack
    tuple val(meta), path("*/expected.copy.number.bed")       , emit: expected_cn               , optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def expected_cn = expected_cn_bed ? "--expected-cn ${expected_cn_bed}" : ""
    def maf = maf_vcf ? "--maf ${maf_vcf}" : ""
    def cnv_exclude_regions = cnv_exclude_regions ? "--cnv-exclude-regions ${cnv_exclude_regions}" : ""

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def expected_cn_bed = expected_cn_bed ? "touch ${prefix}/expected.copy.number.bed": ''
    def maf_mpack = maf_vcf ? "echo \"MAF VCF mpack\" > ${prefix}/maf.mpack": ''
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
    ${expected_cn_bed}
    ${maf_mpack}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sawfish: \$(sawfish --version | sed 's/sawfish //g')
    END_VERSIONS
    """
}
