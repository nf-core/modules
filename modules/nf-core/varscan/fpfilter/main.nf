process VARSCAN_FPFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varscan:2.4.6--hdfd78af_0':
        'biocontainers/varscan:2.4.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(rc)

    output:
    tuple val(meta), path("*.pass.vcf.gz"), emit: pass_vcf
    tuple val(meta), path("*.fail.vcf.gz"), emit: fail_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkfifo vcf_file
    gzip -cdf $vcf > vcf_file &

    varscan fpfilter \\
        vcf_file \\
        $rc \\
        --output-file ${prefix}.varscan.pass.vcf \\
        --filtered-file ${prefix}.varscan.fail.vcf \\
        $args

    gzip ${prefix}.varscan.pass.vcf
    gzip ${prefix}.varscan.fail.vcf

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

    touch ${prefix}.varscan.pass.vcf
    echo "test" > ${prefix}.varscan.pass.vcf
    touch ${prefix}.varscan.fail.vcf
    echo "test" > ${prefix}.varscan.fail.vcf

    gzip ${prefix}.varscan.pass.vcf
    gzip ${prefix}.varscan.fail.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(varscan 2>&1 | grep VarScan | head -n 1 | sed 's/VarScan //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
