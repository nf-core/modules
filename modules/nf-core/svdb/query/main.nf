process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ff995c756aa8a3c0af13b1d054eacd536a11d35de5fa288dacf558bc21696968/data':
        'community.wave.seqera.io/library/bcftools_svdb:ae3b14d2d608fd81' }"

    input:
    tuple val(meta), path(vcf)
    val(in_occs)
    val(in_frqs)
    val(out_occs)
    val(out_frqs)
    path(vcf_dbs)
    path(bedpe_dbs)

    output:
    tuple val(meta), path("*_query.vcf")    , emit: vcf
    tuple val("${task.process}"), val('svdb'), eval("svdb | sed -nE 's/.*SVDB-([0-9.]+).*/\\1/p'"), emit: versions_svdb, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def in_occ           = ""
    def in_frq           = ""
    def out_occ          = ""
    def out_frq          = ""
    def dbs_argument     = vcf_dbs  ? "--db ${vcf_dbs.join(',')}" : ''
    def bedpeds_argument = bedpe_dbs ? "--bedpedb ${bedpe_dbs.join(',')}" : ''
    if (in_occs) {
        in_occ  = "--in_occ ${in_occs.join(',')}"
    }
    if (in_frqs) {
        in_frq  = "--in_frq ${in_frqs.join(',')}"
    }
    if (out_occs) {
        out_occ = "--out_occ ${out_occs.join(',')}"
    }
    if (out_frqs) {
        out_frq = "--out_frq ${out_frqs.join(',')}"
    }
    if ( vcf_dbs && bedpe_dbs ) error "bedpedb input is not compatible with db inputs"
    """
    svdb \\
        --query \\
        $in_occ \\
        $in_frq \\
        $out_occ \\
        $out_frq \\
        $args \\
        $dbs_argument \\
        $bedpeds_argument \\
        --query_vcf $vcf \\
        --prefix ${prefix}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_query.vcf

    """
}
