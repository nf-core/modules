process DEACON_INDEX {
    tag "fasta - ${meta.id} - ${subcommand ?: 'build'}"

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.13.2--h7ef3eeb_1':
        'biocontainers/deacon:0.13.2--h7ef3eeb_0' }"

    input:
    tuple val(meta), path(reference)        // reference fasta file (for build) or idx index (for set operations)
    val(subcommand)                         // build, union, intersect or diff
    tuple val(meta_set), path(set_genomes)  // one (or more for union) idx indexes (or fasta files for diff)

    output:
    tuple val(meta_out), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // default to build if no subcommand is provided
    subcommand = subcommand ?: "build"

    // validate provided subcommand
    def valid_subcommands = ['build', 'union', 'intersect', 'diff']
    if (!(subcommand in valid_subcommands)) {
        error "Invalid deacon index subcommand: '${subcommand}'. Must be one of: ${valid_subcommands.join(', ')}, or empty to default to `build`"
    }

    // create new meta map that stores utilised deacon index subcommand
    meta_out = meta + [operation: subcommand]

    // validate inputs based on subcommand
    if (subcommand == 'build') {
        if (set_genomes) {
            error "Subcommand 'build' does not accept set_genomes input. Received: ${set_genomes}"
        }
    }
    else {
        // union, intersect, or diff require set_genomes input
        if (!set_genomes) {
            error "Subcommand '${subcommand}' requires set_genomes input (idx or fasta files)"
        }

        // create list of files (in case of a single file input)
        def files_to_check = [set_genomes].flatten()

        // union and intersect accept 1 or more idx files
        if (subcommand == 'union' || subcommand == 'intersect') {
            def non_idx = files_to_check.findAll { file -> !file.name.endsWith('.idx') }
            if (non_idx) {
                error "Subcommand '${subcommand}' requires .idx files. Found non-index files: ${non_idx.collect { it.name }.join(', ')}"
            }
        }

        // diff accepts exactly 1 idx or fasta file
        else if (subcommand == 'diff') {
            def valid_exts = ['.idx', '.fa', '.fasta', '.fna']
            def invalid_files = files_to_check.findAll { file ->
                    !valid_exts.any { ext -> file.name.endsWith(ext) }
            }
            if (invalid_files) {
                error "Subcommand 'diff' requires .idx or fasta files (.fa, .fasta, .fna). Found invalid files: ${invalid_files.collect { it.name }.join(', ')}"
            }

            if (files_to_check.size() > 1) {
                error "Subcommand 'diff' expects a single file for set operation. Received ${files_to_check.size()} files."
            }
        }
    }

    """
    deacon \\
        index \\
        ${subcommand} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${subcommand == 'build' ? "${reference}" : "${reference} ${set_genomes}"} > ${prefix}.idx
    """

    stub:
    def prefix = task.ext.prefix ?: "${reference.baseName}"
    subcommand = subcommand ?: 'build'
    meta_out = meta + [operation: subcommand]
    """
    touch ${prefix}.idx
    """
}
