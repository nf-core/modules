process GATK4_DETERMINEGERMLINECONTIGPLOIDY {
    tag "$meta.id"
    label 'process_single'

    if(params.enable_conda){
        error "Conda environments cannot be used for GATK4/DetermineGermlineContigPloidy at the moment. Please use docker or singularity containers."
    }
    container "broadinstitute/gatk:4.3.0.0"

    input:
    tuple val(meta), path(counts), path(bed), path(exclude_beds)
    path(contig_ploidy_table)
    path(ploidy_model)

    output:
    tuple val(meta), path("*-calls.tar.gz") , emit: calls
    tuple val(meta), path("*-model.tar.gz") , emit: model
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_list = counts.collect(){"--input $it"}.join(" ")
    def intervals = bed ? "--intervals ${bed}" : ""
    def exclude = exclude_beds ? exclude_beds.collect(){"--exclude-intervals $it"}.join(" ") : ""
    def model = ploidy_model ? "--model ${ploidy_model}" : ""
    def contig_ploidy = contig_ploidy_table ? "--contig-ploidy-priors ${contig_ploidy_table}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK DetermineGermlineContigPloidy] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" DetermineGermlineContigPloidy \\
        ${input_list} \\
        --output ./ \\
        --output-prefix ${prefix} \\
        ${intervals} \\
        ${exclude} \\
        ${contig_ploidy} \\
        ${model} \\
        --tmp-dir . \\
        ${args}

    tar czf ${prefix}-calls.tar.gz ${prefix}-calls
    tar czf ${prefix}-model.tar.gz ${prefix}-model

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
