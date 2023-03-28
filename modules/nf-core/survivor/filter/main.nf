process SURVIVOR_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::survivor=1.0.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/survivor:1.0.7--h9a82719_1':
        'quay.io/biocontainers/survivor:1.0.7--h9a82719_1' }"

    input:
    tuple val(meta), path(vcf_file) // VCF file to filter
    path(bed)           // file with regions to ignore (NA to disable)
    val(minsv)          // Min SV size (-1 to disable)
    val(maxsv)          // Max SV size (-1 to disable)
    val(minallelefreq)  // Min allele frequency (0-1)
    val(minnumreads)    // Min number of reads support: RE flag (-1 to disable)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_file = bed ? "${bed}" : "NA"

    """
    SURVIVOR \\
        filter \\
        $vcf_file \\
        $bed_file \\
        $minsv \\
        $maxsv \\
        $minallelefreq \\
        $minnumreads \\
        ${prefix}.survivor_filter.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
