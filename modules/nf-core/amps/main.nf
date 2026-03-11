process AMPS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hops:0.35--hdfd78af_1' :
        'biocontainers/hops:0.35--hdfd78af_1' }"

    input:
    path maltextract_results
    path taxon_list
    val filter

    output:
    path "results/heatmap_overview_Wevid.json", emit: json
    path "results/heatmap_overview_Wevid.pdf" , emit: summary_pdf
    path "results/heatmap_overview_Wevid.tsv" , emit: tsv
    path "results/pdf_candidate_profiles/"    , emit: candidate_pdfs
    tuple val("${task.process}"), val('hops'), eval("hops --version 2>&1 | sed 's/HOPS version//' "), emit: versions_hops, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    postprocessing.AMPS.r \\
        -r ${maltextract_results} \\
        -n ${taxon_list} \\
        -m ${filter} \\
        -t ${task.cpus} \\
        -j \\
        ${args}

    """

    stub:

    """
    mkdir -p results/pdf_candidate_profiles
    touch results/heatmap_overview_Wevid.json
    touch results/heatmap_overview_Wevid.pdf
    touch results/heatmap_overview_Wevid.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amps: \$(echo \$(hops --version 2>&1) | sed 's/HOPS version//')
    END_VERSIONS
    """
}
