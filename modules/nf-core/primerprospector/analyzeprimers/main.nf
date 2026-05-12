process PRIMERPROSPECTOR_ANALYZEPRIMERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/primerprospector:1.0.1--py27_0' :
        'quay.io/biocontainers/primerprospector:1.0.1--py27_0' }"

    input:
    tuple val(meta), path(fasta), path(primers)

    output:
    tuple val(meta), path("${prefix}_hits.txt"), emit: hits
    tuple val(meta), path("${prefix}.ps")      , emit: plots
    tuple val("${task.process}"), val('primerprospector'), eval("analyze_primers.py --version 2>&1 | sed 's/.* //; s/-release//'"), topic: versions, emit: versions_primerprospector

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_arg = fasta instanceof List ? fasta.join(':') : fasta
    def primer_arg = primers ? "-P \"${primers}\"" : ''
    def arg_tokens = args.tokenize()
    if (arg_tokens.any { arg -> arg == '-f' || arg.startsWith('--fasta_seqs') || arg == '-P' || arg.startsWith('--primers_filepath') || arg == '-o' || arg.startsWith('--output_dir') }) {
        error "'-f/--fasta_seqs', '-P/--primers_filepath' and '-o/--output_dir' are reserved by this module. Use input files for fasta/primers and task.ext.prefix to control output file prefixes."
    }
    if (!primers && !(arg_tokens.any { arg -> arg == '-p' || arg.startsWith('--primer_name') } && arg_tokens.any { arg -> arg == '-s' || arg.startsWith('--primer_sequence') })) {
        error "Provide a primers file in the input tuple, or specify a single primer with '-p/--primer_name' and '-s/--primer_sequence' in task.ext.args."
    }
    """
    # Primer Prospector 1.0.1 passes numpy floats to range() while plotting.
    # This shim avoids a Python 2 TypeError that otherwise stops .ps output generation.
    cat <<'PY' > analyze_primers_compat.py
    import __builtin__
    _range = __builtin__.range
    def _int_range(*args):
        return _range(*[int(arg) for arg in args])
    __builtin__.range = _int_range
    execfile('/usr/local/bin/analyze_primers.py')
    PY

    python analyze_primers_compat.py \\
        $args \\
        -f "${fasta_arg}" \\
        ${primer_arg}

    mv *_hits.txt "${prefix}_hits.txt"
    mv *.ps "${prefix}.ps"
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def arg_tokens = args.tokenize()
    if (arg_tokens.any { arg -> arg == '-f' || arg.startsWith('--fasta_seqs') || arg == '-P' || arg.startsWith('--primers_filepath') || arg == '-o' || arg.startsWith('--output_dir') }) {
        error "'-f/--fasta_seqs', '-P/--primers_filepath' and '-o/--output_dir' are reserved by this module. Use input files for fasta/primers and task.ext.prefix to control output file prefixes."
    }
    if (!primers && !(arg_tokens.any { arg -> arg == '-p' || arg.startsWith('--primer_name') } && arg_tokens.any { arg -> arg == '-s' || arg.startsWith('--primer_sequence') })) {
        error "Provide a primers file in the input tuple, or specify a single primer with '-p/--primer_name' and '-s/--primer_sequence' in task.ext.args."
    }
    """
    echo "$args"

    touch ${prefix}_hits.txt
    touch ${prefix}.ps
    """
}
