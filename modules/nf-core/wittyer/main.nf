process WITTYER {
    tag "$meta.id"
    label 'process_single'

    container "nf-core/wittyer:0.3.3.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "WITTYER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(query_vcf), path(query_vcf_index), path(truth_vcf), path(truth_vcf_index), path(bed)

    output:
    tuple val(meta),    path("*.json")         , emit: report
    tuple val(meta),    path("*.vcf.gz")       , emit: bench_vcf
    tuple val(meta),    path("*.vcf.gz.tbi")   , emit: bench_vcf_tbi
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--includeBed=$bed" : ""
    if ("$truth_vcf" == "${prefix}.vcf.gz")         error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ("$query_vcf" == "${prefix}.vcf.gz")         error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ("$query_vcf_index" == "${prefix}.vcf.gz.tbi") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ("$query_vcf_index" == "${prefix}.vcf.gz.tbi") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    // dotnet /opt/Wittyer/Wittyer.dll might need to be replaced with new docker image
    """
    mkdir bench

    dotnet /opt/Wittyer/Wittyer.dll \\
        --truthVcf=${truth_vcf} \\
        --inputVcf=${query_vcf} \\
        --outputDirectory=bench \\
        ${regions} \\
        ${args}

    mv bench/Wittyer.Stats.json ${prefix}.json
    mv bench/*.vcf.gz ${prefix}.vcf.gz
    mv bench/*.vcf.gz.tbi ${prefix}.vcf.gz.tbi

    rm -rf bench

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wittyer: \$(dotnet /opt/Wittyer/Wittyer.dll --version  |& sed '1!d ; s/witty.er //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wittyer: \$(dotnet /opt/Wittyer/Wittyer.dll --version  |& sed '1!d ; s/witty.er //')
    END_VERSIONS
    """
}
