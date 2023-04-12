process GATK4_DETERMINEGERMLINECONTIGPLOIDY {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    container "broadinstitute/gatk:4.4.0.0" //Biocontainers is missing a package

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "GATK4_DETERMINEGERMLINECONTIGPLOIDY module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(counts), path(bed), path(exclude_beds)
    path(contig_ploidy_table)
    path(ploidy_model)

    output:
    tuple val(meta), path("*-calls.tar.gz") , emit: calls
    tuple val(meta), path("*-model.tar.gz") , emit: model, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = counts.collect(){"--input $it"}.join(" ")
    def intervals = bed ? "--intervals ${bed}" : ""
    def exclude = exclude_beds ? exclude_beds.collect(){"--exclude-intervals $it"}.join(" ") : ""
    def untar_model = ploidy_model ? (ploidy_model.name.endsWith(".tar.gz") ? "tar -xzf ${ploidy_model}" : "") : ""
    def tar_model = ploidy_model ? "" : "tar czf ${prefix}-model.tar.gz ${prefix}-model"
    def model = ploidy_model ? (ploidy_model.name.endsWith(".tar.gz") ? "--model ${ploidy_model.toString().replace(".tar.gz","")}" : "--model ${ploidy_model}") : ""
    def contig_ploidy = contig_ploidy_table ? "--contig-ploidy-priors ${contig_ploidy_table}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK DetermineGermlineContigPloidy] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    ${untar_model}

    gatk --java-options "-Xmx${avail_mem}M" DetermineGermlineContigPloidy \\
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
    ${tar_model}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-calls.tar.gz
    touch ${prefix}-model.tar.gz
    touch ${prefix}.tsv
    touch ${prefix}2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
