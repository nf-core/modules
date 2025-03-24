process PLINK2_VCF2BGEN {
    tag "${meta.id}"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0'
        : 'biocontainers/plink2:2.00a5.10--h4ac6f70_0'}"

    input:
    tuple val(meta), path(vcf), val(dosage_field), val(bgen_reffirst), val(sample_name_mode)

    output:
    tuple val(meta), path("*.bgen"), emit: bgen_file
    tuple val(meta), path("*.sample"), emit: sample_file
    tuple val(meta), path("*.log"), emit: log_file
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reffirst = bgen_reffirst || task.ext.reffirst ? "ref-first" : ""
    def dosage = task.ext.dosage_field ? task.ext.dosage_field : dosage_field
    def sample_name_opt = task.ext.sample_name_mode ? task.ext.sample_name_mode : sample_name_mode
    """
    plink2 \
        --threads ${task.cpus} \
        --memory ${task.memory.toMega()} \
        --vcf ${vcf} 'dosage=${dosage}' \
        --export bgen-1.2 ${reffirst}\
        --max-alleles 2 \
        --vcf-half-call r \
        --${sample_name_opt} \
        --out ${prefix} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bgen
    touch ${prefix}.sample
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
