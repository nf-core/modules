process NGSCHECKMATE_NCM {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a87c8e024fc7d44064c5c304d3b3bd668a88579b9e069d40b74bcc2458d9dc91/data':
        'community.wave.seqera.io/library/bcftools_ngscheckmate:5d4fe3e82ae99a2b' }"

    input:
    tuple val(meta) , path(files)
    tuple val(meta2), path(snp_bed)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*_corr_matrix.txt"), emit: corr_matrix
    tuple val(meta), path("*_matched.txt")    , emit: matched
    tuple val(meta), path("*_all.txt")        , emit: all
    tuple val(meta), path("*.pdf")            , emit: pdf, optional: true
    tuple val(meta), path("*.vcf")            , emit: vcf, optional: true
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    def unzip = files.any { it.toString().endsWith(".vcf.gz") }
    """
    if $unzip
    then
        for VCFGZ in *.vcf.gz; do
            gunzip -cdf \$VCFGZ > \$( basename \$VCFGZ .gz );
        done
    fi

    NCM_REF="./"${fasta} ncm.py -d . -bed ${snp_bed} -O . -N ${prefix} $args

    if $unzip
    then
        rm -f *.vcf  # clean up decompressed vcfs
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    """
    touch ${prefix}_output_corr_matrix.txt
    touch ${prefix}_matched.txt
    touch ${prefix}_all.txt
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """

}
