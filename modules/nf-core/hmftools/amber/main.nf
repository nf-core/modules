process HMFTOOLS_AMBER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-amber:4.2--hdfd78af_0' :
        'biocontainers/hmftools-amber:4.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_fai
    path loci

    output:
    tuple val(meta), path("${prefix}/"), emit: amber_dir
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tumor_id  = meta.tumor_id  ?: "${meta.id}_tumor"
    def normal_id = meta.normal_id ?: "${meta.id}_normal"
    """
    amber \\
        -tumor ${tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        -reference ${normal_id} \\
        -reference_bam ${normal_bam} \\
        -output_dir ${prefix} \\
        -loci ${loci} \\
        -ref_genome ${fasta} \\
        -threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmftools-amber: \$(amber -version 2>&1 | grep -oP '(?<=AMBER version )[\\d.]+')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${meta.id}.amber.baf.tsv.gz
    touch ${prefix}/${meta.id}.amber.baf.vcf.gz
    touch ${prefix}/${meta.id}.amber.contamination.tsv
    touch ${prefix}/${meta.id}.amber.qc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmftools-amber: \$(amber -version 2>&1 | grep -oP '(?<=AMBER version )[\\d.]+')
    END_VERSIONS
    """
}