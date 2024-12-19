process CELLSNP_MODEA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cellsnp-lite:1.2.3--h6141fd1_2' :
        'biocontainers/cellsnp-lite:1.2.3--h6141fd1_2' }"

    input:
    tuple val(meta), path(bam), path(bai), path(region_vcf), path(barcode)

    output:
    tuple val(meta), path("*.base.vcf.gz") , emit: base
    tuple val(meta), path("*.cells.vcf.gz"), emit: cell          , optional: true
    tuple val(meta), path("*.samples.tsv") , emit: sample
    tuple val(meta), path("*.tag.AD.mtx")  , emit: allele_depth
    tuple val(meta), path("*.tag.DP.mtx")  , emit: depth_coverage
    tuple val(meta), path("*.tag.OTH.mtx") , emit: depth_other
    path 'versions.yml'                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def region_file  = region_vcf      ? "-R $region_vcf"  : ''
    """
    cellsnp-lite -s $bam \\
        -b $barcode \\
        $region_file \\
        -O . \\
        --gzip \\
        --nproc $task.cpus \\
        $args

    mv cellSNP.base.vcf.gz ${prefix}.base.vcf.gz
    if [[ "$args" == *"--genotype"* ]]; then
        mv cellSNP.cells.vcf.gz ${prefix}.cells.vcf.gz
    fi
    mv cellSNP.tag.AD.mtx ${prefix}.tag.AD.mtx
    mv cellSNP.tag.DP.mtx ${prefix}.tag.DP.mtx
    mv cellSNP.tag.OTH.mtx ${prefix}.tag.OTH.mtx
    mv cellSNP.samples.tsv ${prefix}.samples.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp: \$(cellsnp-lite --v | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    echo "" | gzip > ${prefix}.base.vcf.gz
    touch ${prefix}.samples.tsv
    touch ${prefix}.tag.AD.mtx
    touch ${prefix}.tag.DP.mtx
    touch ${prefix}.tag.OTH.mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp: \$(cellsnp-lite --v | awk '{print \$2}')
    END_VERSIONS
    """
}
