process DAMAGEPROFILER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/damageprofiler:1.1--hdfd78af_2' :
        'biocontainers/damageprofiler:1.1--hdfd78af_2' }"

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
    def reference    = fasta       ? "-r $fasta"        : ""
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

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p $prefix

    touch $prefix/3pGtoA_freq.txt
    touch $prefix/3p_freq_misincorporations.txt
    touch $prefix/5pCtoT_freq.txt
    touch $prefix/5p_freq_misincorporations.txt
    touch $prefix/DNA_comp_genome.txt
    touch $prefix/DNA_composition_sample.txt
    touch $prefix/DamagePlot.pdf
    touch $prefix/DamagePlot_five_prime.svg
    touch $prefix/DamagePlot_three_prime.svg
    touch $prefix/DamageProfiler.log
    touch $prefix/Length_plot.pdf
    touch $prefix/Length_plot_combined_data.svg
    touch $prefix/Length_plot_forward_reverse_separated.svg
    touch $prefix/dmgprof.json
    touch $prefix/editDistance.txt
    touch $prefix/edit_distance.pdf
    touch $prefix/edit_distance.svg
    touch $prefix/lgdistribution.txt
    touch $prefix/misincorporation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        damageprofiler: \$(damageprofiler -v | sed 's/^DamageProfiler v//')
    END_VERSIONS
    """
}
