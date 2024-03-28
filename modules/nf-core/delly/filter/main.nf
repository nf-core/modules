process DELLY_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delly:1.2.6--hb7e2ac5_0':
        'biocontainers/delly:1.2.6--hb7e2ac5_0' }"

    input:
    tuple val(meta), path(bcf), path(csi)
    val(mode)
    path(samples)
    val(tumor_samples)
    val(control_samples)

    output:
    tuple val(meta), path("${prefix}.{bcf,vcf.gz}")        , emit: bcf
    tuple val(meta), path("${prefix}.{bcf.csi,vcf.gz.tbi}"), emit: csi
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bcf"
    if ("$bcf" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    def samples_arg = ''
    def sample_list = ''
    if (mode == 'somatic') {
        if (!samples && tumor_samples && control_samples) {
            sample_list = tumor_samples.collect {"${it}\ttumor"}.join('\n')
            sample_list += '\n' + control_samples.collect {"${it}\tcontrol"}.join('\n')
        }
        else {
            error "For somatic mode, either lists of sample ids for tumor and control samples or a sample description file must be provided."
        }
        samples_arg = '-s samples.tsv'
    } else if (mode != 'germline') {
        error "Mode must be one of 'germline' or 'somatic'."
    }

    bcf_output = suffix == "bcf" ? "--outfile ${prefix}.bcf" : ""
    vcf_output = suffix == "vcf" ? "| bgzip ${args2} --threads ${task.cpus} --stdout > ${prefix}.vcf.gz && tabix ${args3} ${prefix}.vcf.gz" : ""

    """
    if [ "$samples_arg" ]; then echo "$sample_list" > samples.tsv; fi

    delly \\
        filter \\
        $args \\
        -f $mode \\
        $bcf_output \\
        $samples_arg \\
        $bcf \\
        $vcf_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bcf"

    def bcf_output = suffix == "bcf" ? "touch ${prefix}.bcf && touch ${prefix}.bcf.csi" : ""
    def vcf_output = suffix == "vcf" ? "touch ${prefix}.vcf.gz && touch ${prefix}.vcf.gz.tbi" : ""

    """
    ${bcf_output}
    ${vcf_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}