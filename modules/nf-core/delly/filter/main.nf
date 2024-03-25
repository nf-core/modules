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

    output:
    tuple val(meta), path("${prefix}.bcf")    , emit: bcf
    tuple val(meta), path("${prefix}.bcf.csi"), emit: csi
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bcf" == "${prefix}.bcf") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    def samples_arg = ''
    def sample_list = ''
    if (mode == 'somatic') {
        if (!samples && meta.control_samples && meta.tumor_samples) {
            sample_list = meta.control_samples.collect {"${it}\tcontrol"}.join('\n')
            sample_list += '\n' + meta.tumor_samples.collect {"${it}\ttumor"}.join('\n')
        }
        else {
            error "For somatic mode, a list of sample assignments must be provided."
        }
        samples_arg = '-s samples.tsv'
    } else if (mode != 'germline') {
        error "Mode must be one of 'germline' or 'somatic'."
    }

    """
    if [ "$samples_arg" ]; then echo "$sample_list" > samples.tsv; fi
    
    delly \\
        filter \\
        $args \\
        -f $mode \\
        -o ${prefix}.bcf \\
        $samples_arg \\
        $bcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bcf
    touch ${prefix}.bcf.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}