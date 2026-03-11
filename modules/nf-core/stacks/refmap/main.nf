process STACKS_REFMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stacks:2.68--h077b44d_0':
        'biocontainers/stacks:2.68--h077b44d_0' }"


    input:
    tuple val(meta), path(bam_path)
    path(popmap)

    output:


    tuple val(meta), path("catalog.calls") , emit: catalog_calls
    tuple val(meta), path("catalog.chrs.tsv") , emit: catalog_chrs
    tuple val(meta), path("catalog.fa.gz")  , emit: catalog_fa
    tuple val(meta), path("gstacks.log") , emit: gstacks_log
    tuple val(meta), path("gstacks.log.distribs"), emit: gstacks_log_distribs
    tuple val(meta), path("populations.haplotypes.tsv"), emit: haplotypes
    tuple val(meta), path("populations.hapstats.tsv"), emit: hapstats
    tuple val(meta), path("populations.sumstats.tsv"), emit: sumstats
    tuple val(meta), path("populations.sumstats_summary.tsv"), emit: sumstats_summary
    tuple val(meta), path("populations.log"), emit: populations_log
    tuple val(meta), path("populations.log.distribs"), emit: populations_log_distribs
    tuple val(meta), path("ref_map.log"), emit: ref_map_log
    tuple val(meta), path("populations.snps.vcf"), emit: vcf , optional: true
    tuple val(meta), path("populations.snps.genepop"), emit: genepop , optional: true
    tuple val(meta), path("populations.structure"), emit: structure , optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "refmap_output"
    """
    ref_map.pl \\
        --samples ./ \\
        --popmap ${popmap} \\
        -o . \\
        -T $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stacks: \$(populations -v)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "refmap_output"

    """

    touch catalog.calls
    touch catalog.chrs.tsv
    echo "" | gzip > catalog.fa.gz
    touch gstacks.log
    touch gstacks.log.distribs
    touch populations.haplotypes.tsv
    touch populations.hapstats.tsv
    touch populations.sumstats.tsv
    touch populations.sumstats_summary.tsv
    touch populations.log
    touch populations.log.distribs
    touch ref_map.log
    touch populations.snps.vcf
    touch populations.snps.genepop
    touch populations.structure

    """
}
