process BISCUIT_EPIREAD {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113 bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0':
        'quay.io/biocontainers/mulled-v2-db16f1c237a26ea9245cf9924f858974ff321d6e:17fa66297f088a1bc7560b7b90dc273bf23f2d8c-0' }"

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
    if ("$options_snp_bed" == "${prefix}.bed.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/.bis.amb//'`

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
}
