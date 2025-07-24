process VARSCAN_PROCESSSOMATIC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varscan:2.4.6--hdfd78af_0':
        'biocontainers/varscan:2.4.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.Germline.vcf.gz"), emit: germline_vcf
    tuple val(meta), path("*.Germline.hc.vcf.gz"), emit: germline_hc_vcf
    tuple val(meta), path("*.Somatic.vcf.gz"), emit: somatic_vcf
    tuple val(meta), path("*.Somatic.hc.vcf.gz"), emit: somatic_hc_vcf
    tuple val(meta), path("*.LOH.vcf.gz"), emit: loh_vcf
    tuple val(meta), path("*.LOH.hc.vcf.gz"), emit: loh_hc_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def vcf_basename = vcf.name.replaceAll(/\.gz$/, '')
    """
    gzip -fd $vcf

    varscan processSomatic \\
        $args \\
        $vcf_basename

    gzip *.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(varscan 2>&1 | grep VarScan | head -n 1 | sed 's/VarScan //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    echo "" | gzip > ${prefix}.Germline.vcf.gz
    echo "" | gzip > ${prefix}.Germline.hc.vcf.gz
    echo "" | gzip > ${prefix}.Somatic.vcf.gz
    echo "" | gzip > ${prefix}.Somatic.hc.vcf.gz
    echo "" | gzip > ${prefix}.LOH.vcf.gz
    echo "" | gzip > ${prefix}.LOH.hc.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(varscan 2>&1 | grep VarScan | head -n 1 | sed 's/VarScan //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
