// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FAIRBIOHEADERS_COMBINE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d07f6d7a9ec156ba40d7b5c2b4f2f1685aea26601fdcec7252f3146fc61fceb0/data'
        : 'community.wave.seqera.io/library/fair-bioheaders:0.8.5--e659fa573f6f4d5b' }"

    input:
    tuple val(meta), path(metadata), path(sequence)

    output:
    tuple val(meta), path("${prefix}.{fhr,fhp}.fasta{,.gz}"), emit: fasta, optional: true
    tuple val(meta), path("${prefix}.fhr.gfa{,.gz}")         , emit: gfa  , optional: true
    tuple val(meta), path("${prefix}.fhgff3.gff3{,.gz}")     , emit: gff3 , optional: true
    tuple val("${task.process}"), val('fair-bioheaders'), eval("bioheaders --version"), topic: versions, emit: versions_fairbioheaders

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def quote = { value -> "'" + value.toString().replace("'", "'\\''") + "'" }
    def type = sequence.name ==~ /(?i).*\.(fa|fasta|fna|fas|faa)(\.b?gz)?/ ? 'fasta'
        : sequence.name ==~ /(?i).*\.gfa(\.b?gz)?/ ? 'gfa'
        : sequence.name ==~ /(?i).*\.gff3?(\.b?gz)?/ ? 'gff3'
        : null
    if (!type) {
        error("FAIRBIOHEADERS_COMBINE: unsupported sequence file extension: ${sequence.name}")
    }
    def gz = sequence.name ==~ /(?i).*\.b?gz/ ? '.gz' : ''
    // The toolkit names its output <stem>.<fhr|fhp|fhgff3>.<type>[.gz] from the header type
    // in the metadata, and writes BGZF for .gz.
    def input = "${prefix}.${type}${gz}"
    def stage = input == sequence.name ? '' : "ln -s -- ${quote.call(sequence)} ${quote.call(input)}"
    """
    ${stage}

    bioheaders \\
        combine \\
        ${args} \\
        -- \\
        ${quote.call(metadata)} \\
        ${quote.call(input)}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def quote = { value -> "'" + value.toString().replace("'", "'\\''") + "'" }
    def outputs = ["${prefix}.fhr.fasta", "${prefix}.fhr.gfa", "${prefix}.fhgff3.gff3"]
    def create = sequence.name ==~ /(?i).*\.b?gz/
        ? outputs.collect { output -> "echo '' | gzip > ${quote.call(output + '.gz')}" }.join('\n')
        : "touch ${outputs.collect { output -> quote.call(output) }.join(' ')}"
    """
    ${create}
    """
}
