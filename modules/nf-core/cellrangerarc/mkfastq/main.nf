process CELLRANGERARC_MKFASTQ {
    tag "mkfastq"
    label 'process_medium'

    // WARNING !! Cell Ranger ARC mkfastq results are not deterministic, so the number of threads used in the process might affect the results.

    container "nf-core/cellranger-arc-mkfastq:2.0.2"

    input:
    tuple val(meta), path(bcl)
    path csv

    output:
    tuple val(meta), path("${prefix}/outs/fastq_path/*.fastq.gz"), emit: fastq
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_MKFASTQ module does not support Conda. Please use docker or singularity containers."
    }

    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_mkfastq"
    """
    cellranger-arc mkfastq --id=${prefix} \\
        --localmem=${task.memory.toGiga()} \\
        --localcores=1 \\
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
    echo | gzip > ${prefix}/outs/fastq_path/Undetermined_S0_L001_I1_001.fastq.gz
    echo | gzip > ${prefix}/outs/fastq_path/Undetermined_S0_L001_R1_001.fastq.gz
    echo | gzip > ${prefix}/outs/fastq_path/Undetermined_S0_L001_R2_001.fastq.gz
    echo | gzip > ${prefix}/outs/fastq_path/Undetermined_S0_L001_R3_001.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangerarc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
