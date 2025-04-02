process SEVERUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/severus:1.4--pyhdfd78af_0':
        'biocontainers/severus:1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(target_input), path(target_index), path(control_input), path(control_index), path(vcf)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("${prefix}/severus.log")                              , emit: log
    tuple val(meta), path("${prefix}/read_qual.txt")                            , emit: read_qual
    tuple val(meta), path("${prefix}/breakpoints_double.csv")                   , emit: breakpoints_double
    tuple val(meta), path("${prefix}/read_alignments")                          , emit: read_alignments                  , optional: true
    tuple val(meta), path("${prefix}/read_ids.csv")                             , emit: read_ids                         , optional: true
    tuple val(meta), path("${prefix}/severus_collaped_dup.bed")                 , emit: collapsed_dup                    , optional: true
    tuple val(meta), path("${prefix}/severus_LOH.bed")                          , emit: loh                              , optional: true
    tuple val(meta), path("${prefix}/all_SVs/severus_all.vcf")                  , emit: all_vcf                          , optional: true
    tuple val(meta), path("${prefix}/all_SVs/breakpoints_clusters_list.tsv")    , emit: all_breakpoints_clusters_list    , optional: true
    tuple val(meta), path("${prefix}/all_SVs/breakpoints_clusters.tsv")         , emit: all_breakpoints_clusters         , optional: true
    tuple val(meta), path("${prefix}/all_SVs/plots/severus_*.html")             , emit: all_plots                        , optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/severus_somatic.vcf")          , emit: somatic_vcf                      , optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/breakpoints_clusters_list.tsv"), emit: somatic_breakpoints_clusters_list, optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/breakpoints_clusters.tsv")     , emit: somatic_breakpoints_clusters     , optional: true
    tuple val(meta), path("${prefix}/somatic_SVs/plots/severus_*.html")         , emit: somatic_plots                    , optional: true
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def control = control_input ? "--control-bam ${control_input}" : ""
    def vntr_bed = bed ? "--vntr-bed ${bed}" : ""
    def phasing_vcf = vcf ? "--phasing-vcf ${vcf}" : ""
    """
    severus \\
        $args \\
        --threads $task.cpus \\
        --target-bam $target_input \\
        $vntr_bed \\
        $control \\
        $phasing_vcf \\
        --out-dir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        severus: \$(severus --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/all_SVs/plots
    mkdir -p ${prefix}/somatic_SVs/plots

    touch ${prefix}/severus_collaped_dup.bed
    touch ${prefix}/severus.log
    touch ${prefix}/severus_LOH.bed
    touch ${prefix}/read_alignments
    touch ${prefix}/read_ids.csv
    touch ${prefix}/read_qual.txt
    touch ${prefix}/breakpoints_double.csv
    touch ${prefix}/all_SVs/severus_all.vcf
    touch ${prefix}/all_SVs/breakpoints_clusters_list.tsv
    touch ${prefix}/all_SVs/breakpoints_clusters.tsv
    touch ${prefix}/all_SVs/plots/severus_0.html
    touch ${prefix}/all_SVs/plots/severus_1.html
    touch ${prefix}/somatic_SVs/severus_somatic.vcf
    touch ${prefix}/somatic_SVs/breakpoints_clusters_list.tsv
    touch ${prefix}/somatic_SVs/breakpoints_clusters.tsv
    touch ${prefix}/somatic_SVs/plots/severus_0.html
    touch ${prefix}/somatic_SVs/plots/severus_1.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        severus: \$(severus --version)
    END_VERSIONS
    """
}
