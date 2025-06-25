process METACACHE_BUILD {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metacache:2.5.0--h077b44d_0':
        'biocontainers/metacache:2.5.0--h077b44d_0' }"

    input:
    tuple val(meta), path(genome_files, stageAs: 'genomes/*')
    path(taxonomy, stageAs: 'taxonomy/*')  // optional. Should be [names.dmp, nodes.dmp], plus optionally merged.dmp
    path(seq2taxid)                        // optional

    output:
    tuple val(meta), path('*.meta'), path('*.cache*'), emit: db
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (task.cpus > 1) { log.warn("'metacache build' cannot be parallelized: ignoring task.cpus > 1") }
    def taxonomy_filenames = taxonomy.collect{ f -> f.fileName.name }.sort()
    assert !taxonomy || taxonomy_filenames == ['names.dmp', 'nodes.dmp'] || taxonomy_filenames == ['merged.dmp', 'names.dmp', 'nodes.dmp']
    def taxonomy_args = taxonomy ? "-taxonomy taxonomy" : ''
    def seq2taxid_args = seq2taxid ? "-taxpostmap '${seq2taxid}'" : ''
    """
    metacache \\
        build \\
        ${prefix} \\
        genomes/ \\
        $taxonomy_args \\
        $seq2taxid_args \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacache: \$(metacache info |& sed -n 's/^MetaCache version \\+\\([0-9.]\\+\\).*\$/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (task.cpus > 1) { log.warn("'metacache build' cannot be parallelized: ignoring task.cpus > 1") }
    def taxonomy_filenames = taxonomy.collect{ f -> f.fileName.name }.sort()
    assert !taxonomy || taxonomy_filenames == ['names.dmp', 'nodes.dmp'] || taxonomy_filenames == ['merged.dmp', 'names.dmp', 'nodes.dmp']
    def n_outputs = (args ==~ /-parts\s/) ? (args.replaceAll(/^.*\s-parts\s+(\S+).*$/, '$1') as Integer) : 1
    assert n_outputs > 0
    """
    touch '${prefix}.meta'
    touch '${prefix}.cache'{0..${n_outputs-1}}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacache: \$(metacache info |& sed -n 's/^MetaCache version \\+\\([0-9.]\\+\\).*\$/\\1/p')
    END_VERSIONS
    """
}
