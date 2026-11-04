// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AMPS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/908fb0fd95516f91026df526e5603b82f5b663be977a70a10131aa6b01dcde5b/data'
        : 'community.wave.seqera.io/library/hops:0.35--ec5d92deb3e899f2'}"

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
