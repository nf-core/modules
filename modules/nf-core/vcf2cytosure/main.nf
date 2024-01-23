process VCF2CYTOSURE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf2cytosure:0.9.1--pyh7cba7a3_1':
        'biocontainers/vcf2cytosure:0.9.1--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(sv_vcf)
    tuple val(meta2), path(coverage_bed)
    tuple val(meta3), path(cns)
    tuple val(meta4), path(snv_vcf)
    path(blacklist_bed)

    output:
    tuple val(meta), path("*.cgh"), emit: cgh
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def coverage = coverage_bed ? "--coverage ${coverage_bed}" : ''
    def cnvkit = cns ? ( coverage_bed ? '' : "--cn ${cns}" ) : ''
    def snv = snv_vcf ? ( coverage_bed ? '' : "--snv ${snv_vcf}" ) : ''
    def blacklist = blacklist_bed ? "--blacklist ${blacklist_bed}" : ''
    def prefix = task.ext.prefix ?: sv_vcf ? "${meta.id}" : "${meta3.id}"

    if ( cns && coverage_bed || snv_vcf && coverage_bed ) error "Coverage_bed input is not compatible with cns and snv"

    """
    vcf2cytosure \\
        --vcf $sv_vcf \\
        --out ${prefix}.cgh \\
        $coverage \\
        $cnvkit \\
        $snv \\
        $blacklist \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2cytosure: \$(echo \$(vcf2cytosure --version 2>&1) | sed 's/^.* cytosure //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def coverage = coverage_bed ? "--coverage ${coverage_bed}" : ''
    def cnvkit = cns ? ( coverage_bed ? '' : "--cn ${cns}" ) : ''
    def snv = snv_vcf ? ( coverage_bed ? '' : "--snv ${snv_vcf}" ) : ''
    def blacklist = blacklist_bed ? "--blacklist ${blacklist_bed}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.cgh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2cytosure: \$(echo \$(vcf2cytosure --version 2>&1) | sed 's/^.* cytosure //' )
    END_VERSIONS
    """
}
