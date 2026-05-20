
process CNVKIT_COVERAGE {
    tag "$meta.id"
    label 'process_medium'
   
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.13--pyhdfd78af_0':
        'quay.io/biocontainers/cnvkit:0.9.13--pyhdfd78af_0' }"

    input:

    tuple val(meta), path(alignment_file), path(interval) 
    // alignement file can be BAM or CRAM
    // in cnvkit version 0.9.13, a bedGraph (.bed.gz) file is also accepted. extract per-base coverage from BAM files once (eg. with bedtools genomecov -gb)
    // interval file can be target or antitarget bed file.

    output:
    
    tuple val(meta), path("*.cnn"), emit: coverage
    
    // version
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), topic: versions, emit: versions_cnvkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" // prefix is used to distinguish between target and antitarget files in output. 
    
    """
    cnvkit.py \\
        coverage \\
            ${alignment_file} \\
            ${interval} \\
            ${args} \\
            --processes ${task.cpus} \\
            --output ${prefix}.cnn
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" 
    
    """
    echo $args
    
    touch ${prefix}.cnn

    """
}
