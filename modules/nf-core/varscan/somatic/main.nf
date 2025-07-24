process VARSCAN_SOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varscan:2.4.6--hdfd78af_0':
        'biocontainers/varscan:2.4.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(normal_mpileup), path(tumour_mpileup)

    output:
    tuple val(meta), path("*.snvs.vcf.gz")  , emit: vcf_snvs
    tuple val(meta), path("*.indels.vcf.gz"), emit: vcf_indels
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkfifo normal_in tumour_in
    gzip -cdf $normal_mpileup > normal_in &
    gzip -cdf $tumour_mpileup > tumour_in &

    varscan somatic \\
        normal_in \\
        tumour_in \\
        --output-snp ${prefix}.snvs.vcf \\
        --output-indel ${prefix}.indels.vcf \\
        $args

    rm normal_in tumour_in

    gzip ${prefix}.snvs.vcf
    gzip ${prefix}.indels.vcf

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

    touch ${prefix}.snvs.vcf
    echo "test" > ${prefix}.snvs.vcf
    touch ${prefix}.indels.vcf
    echo "test" > ${prefix}.indels.vcf

    gzip ${prefix}.snvs.vcf
    gzip ${prefix}.indels.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(varscan 2>&1 | grep VarScan | head -n 1 | sed 's/VarScan //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
