process SNPSIFT_DBNSFP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::snpsift=5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:5.1d--hdfd78af_0' :
        'quay.io/biocontainers/snpsift:5.1d--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    tuple val(meta2), path(database), path(dbs_tbi)// TBI files are optional (use when compressed VCF file)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbsparam = task.ext.dbs ? "-f ${task.ext.dbs}": ''
    """
    SnpSift \\
        dbnsfp \\
        $args \\
        $dbsparam \\
        -db $database \\
        $vcf > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """
}
