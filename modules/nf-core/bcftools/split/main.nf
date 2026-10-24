process BCFTOOLS_SPLIT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/164b5421d87e69f8e9f5add5a937e148e8165f11bfb18159da607f56ecdb1495/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--4e7d12da6416a855'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: split_vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools index -s ${vcf} |\\
    cut -f 1 | \\
    while read C; do \\
        bcftools view -O z -o ${prefix}.\${C}.vcf.gz ${vcf} --regions "\${C}" ;\\
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.chr1.vcf.gz"
    touch "${prefix}.chr2.vcf.gz"
    """
}
