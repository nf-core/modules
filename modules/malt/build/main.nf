process MALT_BUILD {

    label 'process_high'

    conda (params.enable_conda ? "bioconda::malt=0.41" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/malt:0.41--1' :
        'quay.io/biocontainers/malt:0.41--1' }"

    input:
    path fastas
    val seq_type
    val classification_type
    path mapping_file
    val mapping_type
    val mapping_db

    output:
    path "malt_index/"   , emit: index
    path "versions.yml"  , emit: versions
    path "malt-build.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6
    if (!task.memory) {
        log.info '[MALT_BUILD] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    def valid_types = ['gi', 'ref', 'syn']
    def valid_dbs = ['eggnog', 'interpro2go', 'kegg', 'seed', 'taxonomy']

    !valid_types.contains(mapping_type) error "Unrecognised mapping_type value for MALT_BUILD. Options: gi, ref, syn"
    !valid_types.contains(mapping_db) error "Unrecognised mapping database value for MALT_BUILD. Options: eggnog, interpro2go, kegg, seed, taxonomy"

    def classification_type = "${mapping_db}" == "taxonomy" ? "Taxonomy" : mapping_db.capitalize()
    def mapping = "--${mapping_type}2${mapping_db} ${mapping_file}"

    """
    malt-build \\
        -v \\
        --input ${fastas.join(' ')} \\
        -s $seq_type \\
        -d 'malt_index/' \\
        -t $task.cpus \\
        -c $classification_type \\
        $args \\
        $mapping |&tee malt-build.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-build --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
