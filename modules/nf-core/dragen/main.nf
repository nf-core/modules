process DRAGEN {
    tag "$meta.id"
    label 'process_long'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'registry.hub.docker.com/etycksen/dragen4:4.2.4' }" // somehow the PATH is not containing dragen, when I run it locally it is there

    input:
    tuple val(meta), path(fastq)
    path reference

    output:
    tuple val(meta) , path("*.bam")                                 , optional: true , emit: bam
    tuple val(meta) , path("*.bai")                                 , optional: true , emit: bai
    tuple val(meta) , path("*.bam.md5sum")                          , optional: true , emit: bam_md5
    tuple val(meta) , path("${prefix}.vcf.gz")                      , optional: true , emit: vcf_gzip
    tuple val(meta) , path("${prefix}.vcf.gz.tbi")                  , optional: true , emit: vcf_gzip_tbi
    tuple val(meta) , path("${prefix}.vcf.gz.md5sum")               , optional: true , emit: vcf_gzip_md5
    tuple val(meta) , path("${prefix}.hard-filtered.vcf.gz")        , optional: true , emit: hard_filtered_vcf_gzip
    tuple val(meta) , path("${prefix}.hard-filtered.vcf.gz.tbi")    , optional: true , emit: hard_filtered_vcf_gzip_tbi
    tuple val(meta) , path("${prefix}.hard-filtered.vcf.gz.md5sum") , optional: true , emit: hard_filtered_vcf_gzip_md5
    tuple val(meta) , path("*_coverage_metrics.csv")                , optional: true , emit: coverage_metrics
    tuple val(meta) , path("*_overall_mean_cov.csv")                , optional: true , emit: overall_mean_cov
    tuple val(meta) , path("*_contig_mean_cov.csv")                 , optional: true , emit: contig_mean_cov
    tuple val(meta) , path("*.insert-stats.tab")                    , optional: true , emit: insert_stats
    tuple val(meta) , path("*.mapping_metrics.csv")                 , optional: true , emit: mapping_metrics
    tuple val(meta) , path("*.g.vcf.gz")                            , optional: true , emit: gvcf_gzip
    tuple val(meta) , path("*.g.vcf.gz.tbi")                        , optional: true , emit: gvcf_gzip_§tbi
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bin_path = task.ext.bin_path ?: 'dragen'

    def input = ''
    // from fastq
    if (fastq) {
        if (fastq.size() > 2) {
            error "Error: cannot have more than 2 fastq files as input."
        } else {
            input = '-1 ' + fastq.join(' -2 ')
        }
    }

    """
    ${bin_path}_reset

    $bin_path \\
        $args \\
        $input \\
        -n $task.cpus \\
        -r $reference \\
        --output-file-prefix $prefix \\
        --output-directory \$(pwd)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$($bin_path --version 2>&1) | sed 's/^dragen Version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo \$($bin_path --version 2>&1) | sed 's/^dragen Version //')
    END_VERSIONS
    """
}
