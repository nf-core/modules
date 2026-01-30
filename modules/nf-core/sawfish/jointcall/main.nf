process SAWFISH_JOINTCALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca71c93b472a8b9a7701a744a5f123e4474ce9bf1e0d110aae5b84b5134dd74c/data' :
        'community.wave.seqera.io/library/sawfish:2.2.0--430c21f2b465b4f7' }"

    input:
    tuple val(meta), path(sample_dirs)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bams), path(bais)
    tuple val(meta4), path(sample_csv)

    output:
    tuple val(meta), path("*/*_genotyped.sv.vcf.gz")               , emit: vcf
    tuple val(meta), path("*/*_genotyped.sv.vcf.gz.tbi")           , emit: tbi
    tuple val(meta), path("*/contig.alignment.bam")                , emit: bam
    tuple val(meta), path("*/contig.alignment.bam.csi")            , emit: bam_index
    tuple val(meta), path("*/run.stats.json")                      , emit: stats
    tuple val(meta), path("*/samples/*/depth.bw")                  , emit: depth_bw
    tuple val(meta), path("*/samples/*/copynum.bedgraph")          , emit: copynum_bedgraph          , optional: true
    tuple val(meta), path("*/samples/*/gc_bias_corrected_depth.bw"), emit: gc_bias_corrected_depth_bw, optional: true
    tuple val(meta), path("*/samples/*/copynum.summary.json")      , emit: copynum_summary           , optional: true
    tuple val(meta), path("*/sawfish.log")                         , emit: log
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_args = sample_csv ? '' : sample_dirs.sort { dir -> dir.name }.collect { dir -> "--sample ${dir}" }.join(' ')
    def sample_csv_arg = sample_csv ? "--sample-csv ${sample_csv}" : ""
    if (sample_args && sample_csv_arg) {
        throw new IllegalArgumentException("--sample-csv cannot be used together with --sample; choose one input method.")
    }
    """
    sawfish \\
        joint-call \\
        --threads $task.cpus \\
        --ref $fasta \\
        $sample_args \\
        $args \\
        $sample_csv_arg \\
        --output-dir ${prefix}

    # Rename the output files to include prefix
    if [ -f ${prefix}/genotyped.sv.vcf.gz ]; then
        mv ${prefix}/genotyped.sv.vcf.gz ${prefix}/${prefix}_genotyped.sv.vcf.gz
    fi
    if [ -f ${prefix}/genotyped.sv.vcf.gz.tbi ]; then
        mv ${prefix}/genotyped.sv.vcf.gz.tbi ${prefix}/${prefix}_genotyped.sv.vcf.gz.tbi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sawfish: \$(sawfish --version | sed 's/sawfish //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/samples/sample0001_test/
    echo \"\" | gzip > ${prefix}/${prefix}_genotyped.sv.vcf.gz
    touch ${prefix}/${prefix}_genotyped.sv.vcf.gz.tbi
    touch ${prefix}/contig.alignment.bam
    touch ${prefix}/contig.alignment.bam.csi
    touch ${prefix}/run.stats.json
    touch ${prefix}/sawfish.log
    touch ${prefix}/samples/sample0001_test/copynum.bedgraph
    touch ${prefix}/samples/sample0001_test/depth.bw
    touch ${prefix}/samples/sample0001_test/gc_bias_corrected_depth.bw
    touch ${prefix}/samples/sample0001_test/copynum.summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sawfish: \$(sawfish --version | sed 's/sawfish //g')
    END_VERSIONS
    """
}
