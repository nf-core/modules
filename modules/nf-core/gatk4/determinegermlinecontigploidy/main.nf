
process GATK4_DETERMINEGERMLINECONTIGPLOIDY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(counts), path(bed), path(exclude_beds)
    tuple val(meta2), path(ploidy_model)
    path(contig_ploidy_table)

    output:
    tuple val(meta), path("${prefix}-calls"), emit: calls
    tuple val(meta), path("${prefix}-model"), emit: model, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args       ?: ''
    prefix            = task.ext.prefix     ?: "${meta.id}"
    def intervals     = bed                 ? "--intervals ${bed}" : ""
    def exclude       = exclude_beds        ? exclude_beds.collect(){"--exclude-intervals $it"}.join(" ") : ""
    def contig_ploidy = contig_ploidy_table ? "--contig-ploidy-priors ${contig_ploidy_table}" : ""
    def model         = ploidy_model        ? "--model ${ploidy_model}" : ""
    def input_list    = counts.collect(){"--input $it"}.join(" ")

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK DetermineGermlineContigPloidy] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    export THEANO_FLAGS="base_compiledir=\$PWD"
    export PYTENSOR_FLAGS="base_compiledir=\$PWD"
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        DetermineGermlineContigPloidy \\
        ${input_list} \\
        --output ./ \\
        --output-prefix ${prefix} \\
        ${intervals} \\
        ${exclude} \\
        ${contig_ploidy} \\
        ${model} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-calls
    touch ${prefix}-model

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
