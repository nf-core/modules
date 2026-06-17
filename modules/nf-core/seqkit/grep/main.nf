process SEQKIT_GREP {
    tag "${meta.id}"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(sequence)
    path pattern
    val out_ext

    output:
    tuple val(meta), path("*.{fa,fq,fa.gz,fq.gz}"), emit: filter
    tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), emit: versions_seqkit, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def compression_suffix = sequence.getExtension() == "gz" ? ".gz" : ""
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def output_type = out_ext ?: ("${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq")
    def suffix = output_type + compression_suffix
    def pattern_file = pattern ? "-f ${pattern}" : ""
    """
    seqkit \\
        grep \\
        ${args} \\
        --threads ${task.cpus} \\
        ${pattern_file} \\
        ${sequence} \\
        -o ${prefix}.${suffix} \\
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def compression_suffix = sequence.getExtension() == "gz" ? ".gz" : ""
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def output_type = out_ext ?: ("${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq")
    def suffix = output_type + compression_suffix
    """
    echo ${args}

    echo "" | gzip > ${prefix}.${suffix}.gz
    """
}
