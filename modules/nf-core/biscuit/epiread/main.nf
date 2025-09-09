process BISCUIT_EPIREAD {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5b542bbe1f99afd494ef07423ea8b52f2b8a081b85f92db2726c283c78da3cf0/data':
        'community.wave.seqera.io/library/biscuit_samtools:84373c8a97fa63b8' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(snp_bed)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(index)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biscuit_cpus = [(task.cpus * 0.9) as int, 1].max()
    def samtools_cpus = (task.cpus - biscuit_cpus < 1) ? biscuit_cpus : (task.cpus - biscuit_cpus)
    // As of 2/25/22, epiread does not support reading a gzipped SNP BED file.
    // This is a bit hacky but allows the user to supply a gzipped OR uncompressed bed file
    def unzip_snp_bed = snp_bed && (snp_bed.toString() =~ /\.gz$/) ? "bgzip -d ${snp_bed}" : ""
    def unzipped_snp_bed = snp_bed ? snp_bed.toString() - ~/\.gz$/: ""
    // SNP BED input is optional
    def options_snp_bed = snp_bed ? "-B ${unzipped_snp_bed}" : ""
    if ("$options_snp_bed" == "${prefix}.bed.gz") error "Input and output names for biscuit epiread are the same, set prefix in module configuration to disambiguate!"
    """
    ln -sf \$(readlink $fasta) $index/$fasta

    $unzip_snp_bed

    biscuit epiread \\
        -@ $biscuit_cpus \\
        $args \\
        $options_snp_bed \\
        $index/$fasta \\
        $bam \\
        | LC_ALL=C sort -k1,1 -k2,2n \\
        | bgzip -@ $samtools_cpus \\
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
    """
    echo | gzip > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
    END_VERSIONS
    """
}
