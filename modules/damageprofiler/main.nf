process DAMAGEPROFILER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::damageprofiler=1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/damageprofiler:1.1--hdfd78af_2' :
        'quay.io/biocontainers/damageprofiler:1.1--hdfd78af_2' }"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fai
    path specieslist

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def reference    = fasta ? "-r $fasta" : ""
    def species_list = specieslist ? "-sf $specieslist" : ""
    """
    damageprofiler \\
        -i $bam \\
        -o $prefix/ \\
        $args \\
        $reference \\
        $species_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        damageprofiler: \$(damageprofiler -v | sed 's/^DamageProfiler v//')
    END_VERSIONS
    """
}
