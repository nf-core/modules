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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_basename = vcf.name.replaceAll(/\.gz$/, '')
    """
    #vcf_uncompressed=\$(basename "$vcf" .gz)

    gzip -d $vcf

    varscan processSomatic\\
        $args \\
        $vcf_basename \\


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

    touch ${prefix}.Germline.vcf
    echo "test" > ${prefix}.Germline.vcf
    touch ${prefix}.Germline.hc.vcf
    echo "test" > ${prefix}.Germline.hc.vcf
    touch ${prefix}.Somatic.vcf
    echo "test" > ${prefix}.Somatic.vcf
    touch ${prefix}.Somatic.hc.vcf
    echo "test" > ${prefix}.Somatic.hc.vcf
    touch ${prefix}.LOH.vcf
    echo "test" > ${prefix}.LOH.vcf
    touch ${prefix}.LOH.hc.vcf
    echo "test" > ${prefix}.LOH.hc.vcf

    gzip ${prefix}.*.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(varscan 2>&1 | grep VarScan | head -n 1 | sed 's/VarScan //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
