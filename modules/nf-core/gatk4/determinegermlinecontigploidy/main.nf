
process GATK4_DETERMINEGERMLINECONTIGPLOIDY {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    container "nf-core/gatk:4.4.0.0" //Biocontainers is missing a package

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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GATK4_DETERMINEGERMLINECONTIGPLOIDY module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GATK4_DETERMINEGERMLINECONTIGPLOIDY module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
