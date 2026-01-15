process TRAITAR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hmmer_prodigal_pandas_pip_pruned:d2b59b60c4871d7e' :
        'community.wave.seqera.io/library/hmmer_prodigal_pandas_pip_pruned:a83f0296374a52e6' }"

    input:
    tuple val(meta), path(proteins)
    val input_type

    output:
    tuple val(meta), path("*/phenotype_prediction/predictions_majority-vote_combined.txt"), emit: predictions_combined
    tuple val(meta), path("*/phenotype_prediction/predictions_single-votes_combined.txt"), emit: predictions_single_votes
    tuple val(meta), path("*/phenotype_prediction/predictions_flat_*.txt"), optional: true, emit: predictions_flat
    tuple val(meta), path("*/predictions_*.txt"), optional: true, emit: predictions_raw
    tuple val(meta), path("*/annotation/pfam/"), optional: true, emit: pfam_annotation
    tuple val(meta), path("*/gene_prediction/"), optional: true, emit: gene_prediction
    tuple val("${task.process}"), val('traitar'), eval('traitar --version 2>&1 | tail -1'), topic: versions, emit: versions_traitar

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_file = proteins.name.endsWith('.gz') ? proteins.name[0..-4] : proteins.name

    // Validate input_type against allowed values
    if (!['from_genes', 'from_proteins', 'from_nucleotides'].contains(input_type)) {
        error("Invalid input_type: ${input_type}. Must be one of: 'from_genes', 'from_proteins', 'from_nucleotides'")
    }

    """
    mkdir -p input_dir pfam_data

    # Decompress if necessary, otherwise link the file
    if [[ "${proteins}" == *.gz ]]; then
        gunzip -c "${proteins}" > "input_dir/${input_file}"
    else
        ln -s "\$(readlink -f ${proteins})" "input_dir/${input_file}"
    fi

    # Create sample file (tab-separated: filename, sample_name)
    cat > samples.txt <<-EOF
    sample_file_name\tsample_name
    ${input_file}\t${meta.id}
    EOF

    # Download PFAM data for traitar annotation
    traitar pfam pfam_data

    # Run traitar phenotype analysis
    traitar phenotype \\
        pfam_data \\
        input_dir \\
        samples.txt \\
        ${input_type} \\
        ${prefix} \
        -c 1 \\
        --overwrite \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/phenotype_prediction
    mkdir -p ${prefix}/annotation/pfam
    mkdir -p ${prefix}/gene_prediction

    touch ${prefix}/phenotype_prediction/predictions_majority-vote_combined.txt
    touch ${prefix}/phenotype_prediction/predictions_single-votes_combined.txt
    touch ${prefix}/annotation/pfam/${meta.id}_filtered_best.dat
    touch ${prefix}/annotation/pfam/${meta.id}_domtblout.dat
    touch ${prefix}/annotation/pfam/summary.dat
    touch ${prefix}/predictions_raw.txt
    touch ${prefix}/predictions_majority-vote.txt
    touch ${prefix}/predictions_single-votes.txt
    touch ${prefix}/predictions_conservative-vote.txt
    touch ${prefix}/phenotype_prediction/predictions_flat_majority-votes_combined.txt
    touch ${prefix}/phenotype_prediction/predictions_flat_single-votes_combined.txt
    touch ${prefix}/gene_prediction/${meta.id}.faa
    touch ${prefix}/gene_prediction/${meta.id}.gff

    # Get version for stub test
    traitar --version 2>&1 | tail -1 > version.txt
    """
}
