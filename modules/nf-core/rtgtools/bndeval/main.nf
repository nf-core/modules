process RTGTOOLS_BNDEVAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0':
        'biocontainers/rtg-tools:3.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(query_vcf), path(query_vcf_tbi), path(truth_vcf), path(truth_vcf_tbi), path(regions_bed)

    output:
    tuple val(meta), path("*.tp.vcf.gz")                , emit: tp_vcf
    tuple val(meta), path("*.tp.vcf.gz.tbi")            , emit: tp_tbi
    tuple val(meta), path("*.fn.vcf.gz")                , emit: fn_vcf
    tuple val(meta), path("*.fn.vcf.gz.tbi")            , emit: fn_tbi
    tuple val(meta), path("*.fp.vcf.gz")                , emit: fp_vcf
    tuple val(meta), path("*.fp.vcf.gz.tbi")            , emit: fp_tbi
    tuple val(meta), path("*.tp-baseline.vcf.gz")       , emit: baseline_vcf
    tuple val(meta), path("*.tp-baseline.vcf.gz.tbi")   , emit: baseline_tbi
    tuple val(meta), path("*.weighted_roc.tsv.gz")      , emit: weighted_roc
    tuple val(meta), path("*.summary.txt")              , emit: summary
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_regions = regions_bed ? "--bed-regions=${regions_bed}" : ""
    def avail_mem = task.memory.toGiga() + "G"

    """
    rtg RTG_MEM=$avail_mem bndeval \\
        ${args} \\
        --baseline=${truth_vcf} \\
        ${bed_regions} \\
        --calls=${query_vcf} \\
        --output=${prefix}

    cd ${prefix}/
    mv done progress ..
    for f in * ; do mv "\$f" "../${prefix}.\$f" ; done
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo | gzip > ${prefix}.tp.vcf.gz
    touch ${prefix}.tp.vcf.gz.tbi
    echo | gzip > ${prefix}.fn.vcf.gz
    touch ${prefix}.fn.vcf.gz.tbi
    echo | gzip > ${prefix}.fp.vcf.gz
    touch ${prefix}.fp.vcf.gz.tbi
    echo | gzip > ${prefix}.tp-baseline.vcf.gz
    touch ${prefix}.tp-baseline.vcf.gz.tbi
    echo | gzip > ${prefix}.weighted_roc.tsv.gz
    touch ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}
