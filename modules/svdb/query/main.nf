process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::svdb=2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.5.0--py39hcbe4a3b_0':
        'quay.io/biocontainers/svdb:2.5.0--py39hcbe4a3b_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(in_occ), val(in_frq), val(out_occ), val(out_frq)
    path (vcf_db)

    output:
    tuple val(meta), path("*_ann_svdbq.vcf"), emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = ""
    if(in_occ) {
            input += " --in_occ ${in_occ}"
    }
    if(in_frq) {
            input += " --in_frq ${in_frq}"
    }
    if(out_occ) {
            input += " --out_occ ${out_occ}"
    }
    if(out_frq) {
            input += " --out_frq ${out_frq}"
    }
    """
    svdb \\
        --query \\
        $input \\
        $args \\
        --db $vcf_db \\
        --query_vcf $vcf \\
        >${prefix}_ann_svdbq.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """
}
