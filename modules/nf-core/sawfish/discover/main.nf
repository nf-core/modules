process SAWFISH_DISCOVER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7f042970a563400d0486306b8c2c681625e35743e50f06ae7e37b515235bbdd/data' :
        'community.wave.seqera.io/library/sawfish:2.0.4--9ccd13315939ebe4' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(expected_cn_bed), path(maf_vcf), path(cnv_exclude_regions)

    output:
    tuple val(meta), path("versions.yml")           , emit: versions
    tuple val(meta), path("*_discover_dir/assembly.regions.bed")           , emit: assembly_regions
    tuple val(meta), path("*_discover_dir/candidate.sv.bcf")               , emit: candidate_sv_bcf
    tuple val(meta), path("*_discover_dir/candidate.sv.bcf.csi")           , emit: candidate_sv_bcf_csi
    tuple val(meta), path("*_discover_dir/contig.alignment.bam")           , emit: contig_alignment_bam
    tuple val(meta), path("*_discover_dir/contig.alignment.bam.csi")       , emit: contig_alignment_bam_csi
    tuple val(meta), path("*_discover_dir/copynum.bedgraph")               , emit: copynum_bedgraph, optional: true
    tuple val(meta), path("*_discover_dir/copynum.mpack")                  , emit: copynum_mpack, optional: true
    tuple val(meta), path("*_discover_dir/debug.breakpoint_clusters.bed")  , emit: debug_breakpoint_clusters
    tuple val(meta), path("*_discover_dir/debug.cluster.refinement.txt")   , emit: debug_cluster_refinement
    tuple val(meta), path("*_discover_dir/discover.settings.json")         , emit: discover_settings
    tuple val(meta), path("*_discover_dir/genome.gclevels.mpack")          , emit: genome_gclevels
    tuple val(meta), path("*_discover_dir/max.depth.bed")                  , emit: max_depth
    tuple val(meta), path("*_discover_dir/run.stats.json")                 , emit: run_stats
    tuple val(meta), path("*_discover_dir/sample.gcbias.mpack")            , emit: sample_gcbias, optional: true
    tuple val(meta), path("*_discover_dir/sawfish.log")                    , emit: log
    tuple val(meta), path("*_discover_dir/depth.mpack")                    , emit: depth_mpack
    tuple val(meta), path("*_discover_dir/expected.copy.number.bed")       , emit: expected_cn, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--ref ${fasta}" : ""
    def expected_cn = expected_cn_bed ? "--expected-cn ${expected_cn_bed}" : ""
    def maf = maf_vcf ? "--maf ${maf_vcf}" : ""
    def cnv_exclude_regions = cnv_exclude_regions ? "--cnv-exclude-regions ${cnv_exclude_regions}" : ""

    """
    sawfish \\
        discover \\
        --ref $fasta \\
        --bam $bam \\
        --threads $task.cpus \\
        $args \\
        $expected_cn \\
        $cnv_exclude_regions \\
        $maf \\
        --output-dir ${prefix}_discover_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sawfish: \$(sawfish --version | sed 's/sawfish //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def expected_cn_bed = expected_cn_bed ? "touch ${prefix}_discover_dir/expected.copy.number.bed": ''
    def maf_mpack = maf_vcf ? "echo \"MAF VCF mpack\" > ${prefix}_discover_dir/maf.mpack": ''
    def copynum_files = !args.contains('--disable-cnv') ? """
        echo "Copynum bedgraph" > ${prefix}_discover_dir/copynum.bedgraph
        echo "Copynum mpack" > ${prefix}_discover_dir/copynum.mpack
        echo "Sample GC bias mpack" > ${prefix}_discover_dir/sample.gcbias.mpack
        echo "Genome GC levels mpack" > ${prefix}_discover_dir/genome.gclevels.mpack
    """ : ''

    """
    mkdir -p ${prefix}_discover_dir
    echo "##fileformat=VCFv4.2" > ${prefix}_discover_dir/candidate.sv.bcf
    echo "VCF index content" > ${prefix}_discover_dir/candidate.sv.bcf.csi
    echo -e "chr1\t1000\t2000" > ${prefix}_discover_dir/assembly.regions.bed
    echo "BAM content" > ${prefix}_discover_dir/contig.alignment.bam
    echo "BAM index content" > ${prefix}_discover_dir/contig.alignment.bam.csi
    echo -e "chr1\t1000\t2000\tcluster1" > ${prefix}_discover_dir/debug.breakpoint_clusters.bed
    echo "Cluster refinement content" > ${prefix}_discover_dir/debug.cluster.refinement.txt
    echo "{ \"setting\": \"value\" }" > ${prefix}_discover_dir/discover.settings.json
    echo -e "chr1\t1000\t50" > ${prefix}_discover_dir/max.depth.bed
    echo "{ \"stat\": \"value\" }" > ${prefix}_discover_dir/run.stats.json
    echo "Log content" > ${prefix}_discover_dir/sawfish.log
    echo "Depth mpack content" > ${prefix}_discover_dir/depth.mpack
    echo "Copy number bedgraph content" > ${prefix}_discover_dir/copynum.bedgraph
    echo "Copy number mpack content" > ${prefix}_discover_dir/copynum.mpack
    echo "GC bias content" > ${prefix}_discover_dir/sample.gcbias.mpack
    echo "Genome GC levels content" > ${prefix}_discover_dir/genome.gclevels.mpack
    echo "Expected copy number content" > ${prefix}_discover_dir/expected.copy.number.bed

    ${copynum_files}
    ${expected_cn_bed}
    ${maf_mpack}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sawfish: \$(sawfish --version | sed 's/sawfish //g')
    END_VERSIONS
    """
}
