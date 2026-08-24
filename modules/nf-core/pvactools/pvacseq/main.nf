process PVACTOOLS_PVACSEQ {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/griffithlab/pvactools:5.3.1"

    input:
    // Required inputs
    tuple val(meta),val(sample_name),val(hla),path(vcf)
    val algorithms                                 // prediction algorithms: e.g., 'all'

    // Optional inputs
    path iedb                                      // IEDB install directory
    path blastp_path                               // blastp binary path
    path genes_of_interest                         // genes of interest file
    path peptide_fasta                             // peptide FASTA file
    path phased_proximal_variants_vcf              // VCF of phased proximal variants

    output:
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.tsv")                                           , optional: true, emit: mhc_i_intermediate_tsv
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.tsv_*")                                         , optional: true, emit: mhc_i_chunks
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.fasta")                                         , optional: true, emit: mhc_i_fasta
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.net_chop.fa")                                   , optional: true, emit: mhc_i_net_chop_fasta
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.all_epitopes.tsv")                              , optional: true, emit: mhc_i_all_epitopes
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.filtered.tsv")                                  , optional: true, emit: mhc_i_filtered
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.all_epitopes.aggregated.tsv")                   , optional: true, emit: mhc_i_all_epitopes_aggregated
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.all_epitopes.aggregated.tsv.reference_matches") , optional: true, emit: mhc_i_reference_matches
    tuple val(meta), path("${prefix}/MHC_Class_I/${sample_name}.all_epitopes.aggregated.metrics.json")          , optional: true, emit: mhc_i_metrics
    tuple val(meta), path("${prefix}/MHC_Class_I/*.R")                                                          , optional: true, emit: mhc_i_r_files
    tuple val(meta), path("${prefix}/MHC_Class_I/www")                                                          , optional: true, emit: mhc_i_www
    tuple val(meta), path("${prefix}/MHC_Class_I/log")                                                          , optional: true, emit: mhc_i_log

    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.tsv")                                          , optional: true, emit: mhc_ii_intermediate_tsv
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.tsv_*")                                        , optional: true, emit: mhc_ii_chunks
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.fasta")                                        , optional: true, emit: mhc_ii_fasta
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.net_chop.fa")                                  , optional: true, emit: mhc_ii_net_chop_fasta
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.all_epitopes.tsv")                             , optional: true, emit: mhc_ii_all_epitopes
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.filtered.tsv")                                 , optional: true, emit: mhc_ii_filtered
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.all_epitopes.aggregated.tsv")                  , optional: true, emit: mhc_ii_all_epitopes_aggregated
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.all_epitopes.aggregated.tsv.reference_matches"), optional: true, emit: mhc_ii_reference_matches
    tuple val(meta), path("${prefix}/MHC_Class_II/${sample_name}.all_epitopes.aggregated.metrics.json")         , optional: true, emit: mhc_ii_metrics
    tuple val(meta), path("${prefix}/MHC_Class_II/*.R")                                                         , optional: true, emit: mhc_ii_r_files
    tuple val(meta), path("${prefix}/MHC_Class_II/www")                                                         , optional: true, emit: mhc_ii_www
    tuple val(meta), path("${prefix}/MHC_Class_II/log")                                                         , optional: true, emit: mhc_ii_log

    tuple val(meta), path("${prefix}/combined/${sample_name}.fasta")                                            , optional: true, emit: combined_fasta
    tuple val(meta), path("${prefix}/combined/${sample_name}.net_chop.fa")                                      , optional: true, emit: combined_net_chop_fasta
    tuple val(meta), path("${prefix}/combined/${sample_name}.all_epitopes.tsv")                                 , optional: true, emit: combined_all_epitopes
    tuple val(meta), path("${prefix}/combined/${sample_name}.filtered.tsv")                                     , optional: true, emit: combined_filtered
    tuple val(meta), path("${prefix}/combined/${sample_name}.all_epitopes.aggregated.tsv")                      , optional: true, emit: combined_all_epitopes_aggregated
    tuple val(meta), path("${prefix}/combined/${sample_name}.all_epitopes.aggregated.tsv.reference_matches")    , optional: true, emit: combined_reference_matches
    tuple val(meta), path("${prefix}/combined/${sample_name}.all_epitopes.aggregated.metrics.json")             , optional: true, emit: combined_metrics
    tuple val(meta), path("${prefix}/combined/*.R")                                                             , optional: true, emit: combined_r_files
    tuple val(meta), path("${prefix}/combined/www")                                                             , optional: true, emit: combined_www

    path "versions.yml"                                                                                         , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // Validation
    assert sample_name : "Error: sample_name must be provided and not empty"
    assert hla : "Error: hla must be provided and not empty"
    assert algorithms : "Error: algorithms must be provided and not empty"
    assert vcf : "Error: vcf file '${vcf}' does not exist or is not a valid path"

    def iedb_opt = iedb ? "--iedb-install-directory ${iedb}" : ""
    def blastp_opt = blastp_path ? "--blastp-path ${blastp_path}" : ""
    def genes_opt = genes_of_interest ? "--genes-of-interest-file ${genes_of_interest}" : ""
    def phased_opt = phased_proximal_variants_vcf ? "-p ${phased_proximal_variants_vcf}" : ""

    // Combine all optional arguments into one string
    def optional_args = [
        iedb_opt,
        blastp_opt,
        genes_opt,
        phased_opt,
        args,
        "-t $task.cpus"
    ].findAll { it }
    .join(" ")
    """
    pvacseq run \\
        $vcf \\
        $sample_name \\
        $hla \\
        $algorithms \\
        $prefix \\
        $optional_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pvactools -v)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = sample_name

    """
    mkdir -p ${prefix}/MHC_Class_I
    mkdir -p ${prefix}/MHC_Class_II
    mkdir -p ${prefix}/combined

    # MHC Class I outputs
    touch ${prefix}/MHC_Class_I/${sample_name}.tsv
    touch ${prefix}/MHC_Class_I/${sample_name}.tsv_1
    touch ${prefix}/MHC_Class_I/${sample_name}.fasta
    touch ${prefix}/MHC_Class_I/${sample_name}.net_chop.fa
    touch ${prefix}/MHC_Class_I/${sample_name}.all_epitopes.tsv
    touch ${prefix}/MHC_Class_I/${sample_name}.filtered.tsv
    touch ${prefix}/MHC_Class_I/${sample_name}.all_epitopes.aggregated.tsv
    touch ${prefix}/MHC_Class_I/${sample_name}.all_epitopes.aggregated.tsv.reference_matches
    touch ${prefix}/MHC_Class_I/${sample_name}.all_epitopes.aggregated.metrics.json
    touch ${prefix}/MHC_Class_I/ui.R
    touch ${prefix}/MHC_Class_I/app.R
    touch ${prefix}/MHC_Class_I/server.R
    touch ${prefix}/MHC_Class_I/styling.R
    touch ${prefix}/MHC_Class_I/anchor_and_helper_functions.R
    mkdir -p ${prefix}/MHC_Class_I/www
    touch ${prefix}/MHC_Class_I/www/anchor.jpg
    touch ${prefix}/MHC_Class_I/www/pVACciew_logo.png
    touch ${prefix}/MHC_Class_I/www/pVACview_logo_mini.png
    mkdir -p ${prefix}/MHC_Class_I/log
    touch ${prefix}/MHC_Class_I/log/inputs.yml

    # MHC Class II outputs
    touch ${prefix}/MHC_Class_II/${sample_name}.tsv
    touch ${prefix}/MHC_Class_II/${sample_name}.tsv_1
    touch ${prefix}/MHC_Class_II/${sample_name}.fasta
    touch ${prefix}/MHC_Class_II/${sample_name}.net_chop.fa
    touch ${prefix}/MHC_Class_II/${sample_name}.all_epitopes.tsv
    touch ${prefix}/MHC_Class_II/${sample_name}.filtered.tsv
    touch ${prefix}/MHC_Class_II/${sample_name}.all_epitopes.aggregated.tsv
    touch ${prefix}/MHC_Class_II/${sample_name}.all_epitopes.aggregated.tsv.reference_matches
    touch ${prefix}/MHC_Class_II/${sample_name}.all_epitopes.aggregated.metrics.json
    touch ${prefix}/MHC_Class_II/ui.R
    touch ${prefix}/MHC_Class_II/app.R
    touch ${prefix}/MHC_Class_II/server.R
    touch ${prefix}/MHC_Class_II/styling.R
    touch ${prefix}/MHC_Class_II/anchor_and_helper_functions.R
    mkdir -p ${prefix}/MHC_Class_II/www
    touch ${prefix}/MHC_Class_II/www/anchor.jpg
    touch ${prefix}/MHC_Class_II/www/pVACciew_logo.png
    touch ${prefix}/MHC_Class_II/www/pVACview_logo_mini.png
    mkdir -p ${prefix}/MHC_Class_II/log
    touch ${prefix}/MHC_Class_II/log/inputs.yml

    # Combined outputs
    touch ${prefix}/combined/${sample_name}.fasta
    touch ${prefix}/combined/${sample_name}.net_chop.fa
    touch ${prefix}/combined/${sample_name}.all_epitopes.tsv
    touch ${prefix}/combined/${sample_name}.filtered.tsv
    touch ${prefix}/combined/${sample_name}.all_epitopes.aggregated.tsv
    touch ${prefix}/combined/${sample_name}.all_epitopes.aggregated.tsv.reference_matches
    touch ${prefix}/combined/${sample_name}.all_epitopes.aggregated.metrics.json
    touch ${prefix}/combined/ui.R
    touch ${prefix}/combined/app.R
    touch ${prefix}/combined/server.R
    touch ${prefix}/combined/styling.R
    touch ${prefix}/combined/anchor_and_helper_functions.R
    mkdir -p ${prefix}/combined/www
    touch ${prefix}/combined/www/anchor.jpg
    touch ${prefix}/combined/www/pVACciew_logo.png
    touch ${prefix}/combined/www/pVACview_logo_mini.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pvactools -v)
    END_VERSIONS
    """
}
