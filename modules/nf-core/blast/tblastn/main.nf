process BLAST_TBLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml"           , emit: versions

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

    DB=`find -L ./ -name "*.nal" | sed 's/\\.nal\$//'`
    if [ -z "\$DB" ]; then
        DB=`find -L ./ -name "*.nin" | sed 's/\\.nin\$//'`
    fi
    echo Using \$DB

    tblastn \\
        -num_threads ${task.cpus} \\
        -db \$DB \\
        -query ${fasta_name} \\
        ${args} \\
        -out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
