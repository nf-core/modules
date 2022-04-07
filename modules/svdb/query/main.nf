process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::svdb=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.6.0--py39h5371cbf_0':
        'quay.io/biocontainers/svdb:2.6.0--py39h5371cbf_0' }"

    input:
    tuple val(meta), path(vcf)
    val(in_occs)
    val(in_frqs)
    val(out_occs)
    val(out_frqs)
    path (vcf_dbs)

    output:
    tuple val(meta), path("*_query.vcf"), emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def in_occ  = ""
    def in_frq  = ""
    def out_occ = ""
    def out_frq = ""
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

    """
    svdb \\
        --query \\
        $in_occ \\
        $in_frq \\
        $out_occ \\
        $out_frq \\
        $args \\
        --db ${vcf_dbs.join(',')} \\
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
