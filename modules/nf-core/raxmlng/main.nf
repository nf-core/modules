process RAXMLNG {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.2.2--h6747034_1' :
        'quay.io/biocontainers/raxml-ng:1.2.2--h6747034_1' }"

    input:
    tuple val(meta), path(alignment), val(model)

    output:
    // either bestTree or bootstraps file is created, depending on options given
    tuple val(meta), path("*.raxml.bestTree")  , emit: phylogeny             , optional:true
    tuple val(meta), path("*.raxml.bootstraps"), emit: phylogeny_bootstrapped, optional:true
    tuple val("${task.process}"), val('raxmlng'), eval("raxml-ng --version 2>&1 | sed '/RAxML-NG v/!d;s/.*v. //;s/ .*//'"), emit: versions_raxmlng, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fix random seed for reproducibility if not specified in command line
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=42"}
    """
    raxml-ng \\
        ${args} \\
        --msa ${alignment} \\
        --model ${model} \\
        --threads ${task.cpus} \\
        --prefix ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_files = args.contains('--bootstrap') || args.contains('--bs-trees') ? "touch ${prefix}.raxml.bootstraps" : "touch ${prefix}.raxml.bestTree"
    """
    ${touch_files}
    """
}
