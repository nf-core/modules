process MITOHIFI_MITOHIFI {
    tag "$meta.id"
    label 'process_high'

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    tuple val(meta) , path(input)
    tuple val(meta2), path(ref_fa)
    tuple val(meta3), path(ref_gb)
    val input_mode
    val mito_code

    output:
    tuple val(meta), path("final_mitogenome.fasta")         , emit: fasta
    tuple val(meta), path("contigs_stats.tsv")              , emit: stats
    tuple val(meta), path("final_mitogenome.gb")            , emit: gb                         , optional: true
    tuple val(meta), path("final_mitogenome.gff")           , emit: gff                        , optional: true
    tuple val(meta), path("all_potential_contigs.fa")       , emit: all_potential_contigs      , optional: true
    tuple val(meta), path("contigs_annotations.png")        , emit: contigs_annotations        , optional: true
    tuple val(meta), path("contigs_circularization/")       , emit: contigs_circularization    , optional: true
    tuple val(meta), path("contigs_filtering/")             , emit: contigs_filtering          , optional: true
    tuple val(meta), path("coverage_mapping/")              , emit: coverage_mapping           , optional: true
    tuple val(meta), path("coverage_plot.png")              , emit: coverage_plot              , optional: true
    tuple val(meta), path("final_mitogenome.annotation.png"), emit: final_mitogenome_annotation, optional: true
    tuple val(meta), path("final_mitogenome_choice/")       , emit: final_mitogenome_choice    , optional: true
    tuple val(meta), path("final_mitogenome.coverage.png")  , emit: final_mitogenome_coverage  , optional: true
    tuple val(meta), path("potential_contigs/")             , emit: potential_contigs          , optional: true
    tuple val(meta), path("reads_mapping_and_assembly/")    , emit: reads_mapping_and_assembly , optional: true
    tuple val(meta), path("shared_genes.tsv")               , emit: shared_genes               , optional: true
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }

    def args      = task.ext.args ?: ''
    def input_arg = ""
    if(input_mode == "contigs") {
        input_arg = "-c ${input}"
    } else if(input_mode == "reads"){
        input_arg = "-r ${input}"
    } else {
        error "Error: invalid MitoHiFi input mode: ${input_mode}. Must be either 'contigs' or 'reads'!"
    }
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION   = '3.2.3'
    """
    mitohifi.py \\
        ${input_arg} \\
        -f ${ref_fa} \\
        -g ${ref_gb} \\
        -o ${mito_code} \\
        -t $task.cpus \\
        ${args}

    ## old version command: \$(mitohifi.py -v | sed 's/.* //')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: ${VERSION}
    END_VERSIONS
    """

    stub:
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '3.2.3'
    """
    touch final_mitogenome.fasta
    touch final_mitogenome.gb
    touch final_mitogenome.gff
    touch contigs_stats.tsv
    touch all_potential_contigs.fa
    touch contigs_annotations.png
    touch coverage_plot.png
    touch final_mitogenome.annotation.png
    touch final_mitogenome.coverage.png
    touch shared_genes.tsv

    mkdir contigs_circularization
    mkdir contigs_filtering
    mkdir coverage_mapping
    mkdir final_mitogenome_choice
    mkdir potential_contigs
    mkdir reads_mapping_and_assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: ${VERSION}
    END_VERSIONS
    """
}
