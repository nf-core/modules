process DIAMOND_BLASTX {
    tag "$meta.id"
    label 'process_medium'

    // Dimaond is limited to v2.0.9 because there is not a
    // singularity version higher than this at the current time.
    conda (params.enable_conda ? "bioconda::diamond=2.0.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0' :
        'quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    val outext

    output:
    tuple val(meta), path('*.{blast,xml,txt,daa,sam,tsv,paf}'), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    switch ( outext ) {
        case "blast": outfmt = 0; break
        case "xml": outfmt = 5; break
        case "txt": outfmt = 6; break
        case "daa": outfmt = 100; break
        case "sam": outfmt = 101; break
        case "tsv": outfmt = 102; break
        case "paf": outfmt = 103; break
    }
    """
    DB=`find -L ./ -name "*.dmnd" | sed 's/.dmnd//'`

    diamond \\
        blastx \\
        --threads $task.cpus \\
        --db \$DB \\
        --query $fasta \\
        --outfmt ${outfmt} \\
        $args \\
        --out ${prefix}.${outext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
