process CELLRANGERARC_MKFASTQ {
    tag "mkfastq"
    label 'process_medium'

    container "nf-core/cellranger-arc-mkfastq:2.0.2"

    input:
    tuple val(meta), path(bcl)
    path csv

    output:
    tuple val(meta), path("versions.yml")                        , emit: versions
    tuple val(meta), path("${prefix}/outs/fastq_path/*.fastq.gz"), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_MKFASTQ module does not support Conda. Please use docker or singularity containers."
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_mkfastq"
    """
    cellranger-arc mkfastq --id=${prefix} \\
        --localmem=${task.memory.toGiga()} \\
        --localcores=${task.cpus} \\
        --run=${bcl} \\
        --csv=${csv} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangerarc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_MKFASTQ module does not support Conda. Please use docker or singularity containers."
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "${prefix}/outs/fastq_path/"
    echo | gzip > ${prefix}/outs/fastq_path/fake_file.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangerarc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
