process BLAST_TBLASTX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5222a42b366a0468a4c795f5057c2b8cfe39489548f8bd807e8ac0f80069bad5/data':
        'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path('*.tblastx.txt'), emit: txt
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`
    if [ -z "\$DB" ]; then
        DB= `find -L ./ -name "*.nin" | sed 's/\\.nin\$//'`
    fi
    echo Using \$DB

    tblastx \\
        -num_threads ${task.cpus} \\
        -db \$DB \\
        -query ${fasta_name} \\
        $args \\
        -out ${prefix}.tblastx.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastx -version 2>&1 | sed 's/^.*tblastx: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tblastx.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastx -version 2>&1 | sed 's/^.*tblastx: //; s/ .*\$//')
    END_VERSIONS
    """
}
