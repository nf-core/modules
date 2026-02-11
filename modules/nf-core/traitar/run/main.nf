process TRAITAR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/hmmer_prodigal_pandas_parallel_pruned:ccae2eabc2a54ac8' :
        'community.wave.seqera.io/library/hmmer_prodigal_pandas_parallel_pruned:ccae2eabc2a54ac8' }"

    input:
    tuple val(meta), path(fasta)
    val input_type
    path pfam_db

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
    def input_file = fasta.name.endsWith('.gz') ? fasta.name[0..-4] : fasta.name

    // Validate input_type against allowed values
    if (!['from_genes', 'from_nucleotides', 'from_annotation_summary'].contains(input_type)) {
        error("Invalid input_type: ${input_type}. Must be one of: 'from_genes', 'from_nucleotides', 'from_annotation_summary'")
    }

    """
    mkdir -p input_dir

    # Decompress if necessary, otherwise link the file
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c "${fasta}" > "input_dir/${input_file}"
    else
        mkdir -p input_dir
        ln -s "\$(readlink -f ${fasta})" "input_dir/${input_file}"
    fi

    # Create sample file (tab-separated: filename, sample_name)
    cat > samples.txt <<-EOF
    sample_file_name\tsample_name
    ${input_file}\t${prefix}
    EOF


    # Run traitar phenotype analysis
    traitar phenotype \\
        ${pfam_db} \\
        input_dir \\
        samples.txt \\
        ${input_type} \\
        ${prefix} \\
        -c ${task.cpus} \\
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
    touch ${prefix}/annotation/pfam/${prefix}_filtered_best.dat
    touch ${prefix}/annotation/pfam/${prefix}_domtblout.dat
    touch ${prefix}/annotation/pfam/summary.dat
    touch ${prefix}/predictions_raw.txt
    touch ${prefix}/predictions_majority-vote.txt
    touch ${prefix}/predictions_single-votes.txt
    touch ${prefix}/predictions_conservative-vote.txt
    touch ${prefix}/phenotype_prediction/predictions_flat_majority-votes_combined.txt
    touch ${prefix}/phenotype_prediction/predictions_flat_single-votes_combined.txt
    # Stub gene prediction with headers
    echo -e ">gene_001\\nMETLQKSTVVA" > ${prefix}/gene_prediction/${prefix}.faa
    echo -e "##gff-version 3\\ngene_001\\t.\\tCDS\\t1\\t30\\t.\\t+\\t0\\t." > ${prefix}/gene_prediction/${prefix}.gff

    # Stub PFAM annotation files with headers
    echo "target_name accession tlen query_name accession qlen evalue bitscore" > ${prefix}/annotation/pfam/${prefix}_domtblout.dat
    echo "accession description" > ${prefix}/annotation/pfam/${prefix}_filtered_best.dat
    echo "done" > ${prefix}/annotation/pfam/summary.dat

    # Add minimal content to phenotype and gene prediction files
    echo -e "header1\\theader2\\n${prefix}\\tvalue" > ${prefix}/phenotype_prediction/predictions_majority-vote_combined.txt
    echo -e "header1\\theader2\\n${prefix}\\tvalue" > ${prefix}/phenotype_prediction/predictions_single-votes_combined.txt
    echo -e ">gene1\\nMKTL" > ${prefix}/gene_prediction/${prefix}.faa
    echo -e "domain1\\tdescription" > ${prefix}/annotation/pfam/${prefix}_filtered_best.dat

    """
}
