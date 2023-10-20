process EGGNOGMAPPER {
    tag "$meta.id"
    label 'process_long'

    conda "bioconda::eggnog-mapper=2.1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(eggnog_db)
    path(eggnog_data_dir)
    path(eggnog_diamond_db)

    output:
    tuple val(meta), path("*.emapper.annotations")   , emit: annotations
    tuple val(meta), path("*.emapper.seed_orthologs"), emit: orthologs
    tuple val(meta), path("*.emapper.hits")          , emit: hits
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")
    def dbmem = task.memory.toMega() > 40000 ? '--dbmem' : ''
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    emapper.py \\
        --cpu ${task.cpus} \\
        -i ${fasta_name} \\
        --data_dir ${eggnog_data_dir} \\
        -m diamond \\
        --dmnd_db ${eggnog_diamond_db} \\
        --database ${eggnog_db} \\
        --output ${prefix} \\
        ${dbmem} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.emapper.annotations
    touch ${prefix}.emapper.seed_orthologs
    touch ${prefix}.emapper.hits

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
    END_VERSIONS
    """
}
