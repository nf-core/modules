// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process METAPHLAN_METAPHLAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1a05e0cb9a7387c4711a0bed4e283152bba02d06712d97c0253aa481b48e1a3e/data'
        : 'community.wave.seqera.io/library/metaphlan:4.3--6cd977ac2f320230'}"

    input:
    tuple val(meta), path(input, stageAs: 'input/*')
    tuple val(meta2), path(db_metaphlan_latest)
    tuple val(meta3), path(ignore_markers)
    val save_samfile

    output:
    tuple val(meta), path("${prefix}_profile.txt"), optional: true, emit: profile
    tuple val(meta), path("${prefix}.biom"), optional: true, emit: biom
    tuple val(meta), path("${prefix}.mapout.bz2"), optional: true, emit: mapout
    tuple val(meta), path("${prefix}.sam.bz2"), optional: true, emit: sam
    tuple val(meta), path("${prefix}_vsc.tsv"), optional: true, emit: vsc
    tuple val(meta), path("${prefix}.subsampled.${input_type}.gz"), optional: true, emit: subsampled
    tuple val("${task.process}"), val('metaphlan'), eval("metaphlan --version 2>&1 | head -n 1 | cut -d ' ' -f 3"), emit: versions_metaphlan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    input_type = "${input}" =~ /.*\.(fastq|fq)/ ? "fastq" : "${input}" =~ /.*\.(fasta|fna|fa)/ ? "fasta" : "${input}" =~ /.*\.mapout(\.bz2)?$/ ? "mapout" : "sam"
    def input_data = input_type == "fastq" && !meta.single_end ? "${input[0]},${input[1]}" : "${input}"
    def is_reads = input_type in ["fastq", "fasta"]
    def profile_vsc = args.contains('--profile_vsc')
    def subsampling = (args =~ /(^|\s)--subsampling(\s|=)/).find() && !args.contains('--mapping_subsampling')
    def map_out = is_reads ? "--mapout ${prefix}.mapout.bz2" : ''
    def samfile_out = (save_samfile || profile_vsc) && is_reads ? "-s ${prefix}.sam.bz2" : ''
    def vsc_out = profile_vsc ? "--vsc_out ${prefix}_vsc.tsv" : ''
    def subsampling_out = subsampling && is_reads ? "--subsampling_output ${prefix}.subsampled.${input_type}.gz" : ''
    def ignore_markers_arg = ignore_markers ? "--ignore_markers ${ignore_markers}" : ''
    def index = args =~ /(^|\s)(--index|-x)\s/ ? '' : "--index \$BT2_DB_INDEX"
    def output_file = args.contains('--biom_format_output') ? "${prefix}.biom" : "${prefix}_profile.txt"
    """
    BT2_DB=`find -L "${db_metaphlan_latest}" -name "*rev.1.bt2*" -exec dirname {} \\; | head -n 1`
    BT2_DB_INDEX=`find -L "${db_metaphlan_latest}" -name "*.rev.1.bt2*" | sed 's/\\.rev.1.bt2.*\$//' | sed 's/.*\\///' | head -n 1`

    metaphlan \\
        --nproc ${task.cpus} \\
        --tmp_dir . \\
        --input_type ${input_type} \\
        ${input_data} \\
        ${args} \\
        ${map_out} \\
        ${samfile_out} \\
        ${vsc_out} \\
        ${subsampling_out} \\
        ${ignore_markers_arg} \\
        --db_dir \$BT2_DB \\
        ${index} \\
        --output_file ${output_file}

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    input_type = "${input}" =~ /.*\.(fastq|fq)/ ? "fastq" : "${input}" =~ /.*\.(fasta|fna|fa)/ ? "fasta" : "${input}" =~ /.*\.mapout(\.bz2)?$/ ? "mapout" : "sam"
    def is_reads = input_type in ["fastq", "fasta"]
    def profile_vsc = args.contains('--profile_vsc')
    def subsampling = (args =~ /(^|\s)--subsampling(\s|=)/).find() && !args.contains('--mapping_subsampling')
    def output_file = args.contains('--biom_format_output') ? "${prefix}.biom" : "${prefix}_profile.txt"
    def map_cmd = is_reads ? "echo '' | bzip2 > ${prefix}.mapout.bz2" : ''
    def samfile_cmd = (save_samfile || profile_vsc) && is_reads ? "echo '' | bzip2 > ${prefix}.sam.bz2" : ''
    def vsc_cmd = profile_vsc ? "touch ${prefix}_vsc.tsv" : ''
    def subsampling_cmd = subsampling && is_reads ? "echo '' | gzip > ${prefix}.subsampled.${input_type}.gz" : ''
    """
    echo "${args}"
    touch ${output_file}
    ${map_cmd}
    ${samfile_cmd}
    ${vsc_cmd}
    ${subsampling_cmd}

    """
}
