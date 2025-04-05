process RAXMLNG {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.0.3--h32fcf60_0' :
        'biocontainers/raxml-ng:1.0.3--h32fcf60_0' }"

    input:
    tuple val(meta), path(alignment), val(model)

    output:
    tuple val(meta), path("*.raxml.bestTree")              , emit: phylogeny
    tuple val(meta), path("*.raxml.support"), optional:true, emit: phylogeny_bootstrapped
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    raxml-ng \\
        $args \\
        --msa $alignment \\
        --model $model \\
        --threads $task.cpus \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: params.raxmlng_args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    if (!meta.id) {
        error "Input meta map does not contain 'id'. Received: ${meta}"
    }
    // Use a dedicated param for stub testing the bootstrap output scenario
    def touch_support = args.contains('--bootstrap') || args.contains('--bs-trees') ? "touch ${prefix}.raxml.support" : ""
    """
    # Create stub output files
    touch ${prefix}.raxml.bestTree
    ${touch_support}

    # Create versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: stub_version
    END_VERSIONS
    """
}
