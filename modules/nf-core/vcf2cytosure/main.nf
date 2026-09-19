process VCF2CYTOSURE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4dc6c9fbe153af0ea2ccb2ff21592a4bf766a676200c6f7efca67cca8e36419f/data':
        'community.wave.seqera.io/library/pip_vcf2cytosure:19948ad2c2df9484' }"

    input:
    tuple val(meta), path(sv_vcf)
    tuple val(meta2), path(coverage_bed)
    tuple val(meta3), path(cns)
    tuple val(meta4), path(snv_vcf)
    path(blacklist_bed)

    output:
    tuple val(meta), path("*.cgh"), emit: cgh
    tuple val("${task.process}"), val('vcf2cytosure'), eval("vcf2cytosure --version 2>&1 | sed -n 's/.*cytosure //p'"), topic: versions, emit: versions_vcf2cytosure

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cgh
    """
}
