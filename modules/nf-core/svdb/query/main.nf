process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.8.2--py39h5371cbf_0':
        'biocontainers/svdb:2.8.2--py39h5371cbf_0' }"

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
    path "versions.yml"                     , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_query.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """
}
