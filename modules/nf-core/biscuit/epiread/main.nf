process BISCUIT_EPIREAD {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0':
        'biocontainers/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(snp_bed)
    path(index)

    output:
    tuple val(meta), path("*.bed.gz"), emit: epiread_bed
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biscuit_cpus = (int) Math.max(Math.floor(task.cpus*0.9),1)
    def samtools_cpus = task.cpus-biscuit_cpus
    // As of 2/25/22, epiread does not support reading a gzipped SNP BED file.
    // This is a bit hacky but allows the user to supply a gzipped OR uncompressed bed file
    def unzip_snp_bed = snp_bed && (snp_bed.toString() =~ /\.gz$/) ? "bgzip -d ${snp_bed}" : ""
    def unzipped_snp_bed = snp_bed ? snp_bed.toString() - ~/\.gz$/: ""
    // SNP BED input is optional
    def options_snp_bed = snp_bed ? "-B ${unzipped_snp_bed}" : ""
    if ("$options_snp_bed" == "${prefix}.bed.gz") error "Input and output names for biscuit epiread are the same, set prefix in module configuration to disambiguate!"
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/\\.bis.amb\$//'`

    $unzip_snp_bed

    biscuit epiread \\
        -@ $biscuit_cpus \\
        $args \\
        $options_snp_bed \\
        \$INDEX \\
        $bam | \\
    LC_ALL=C sort -k1,1 -k2,2n | \\
    bgzip \\
        -@ $samtools_cpus \\
        $args2 \\
        -c > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unzipped_snp_bed = snp_bed ? snp_bed.toString() - ~/\.gz$/: ""
    def options_snp_bed = snp_bed ? "-B ${unzipped_snp_bed}" : ""
    if ("$options_snp_bed" == "${prefix}.bed.gz") error "Input and output names for biscuit epiread are the same, set prefix in module configuration to disambiguate!"
    """
    echo | gzip > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
    END_VERSIONS
    """
}
