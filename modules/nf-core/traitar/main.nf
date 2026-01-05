process TRAITAR {
    tag "$meta.id"
    label 'process_medium'

    container null

    input:
    tuple val(meta), path(proteins)

    output:
    tuple val(meta), path("traitar_out/phenotype_prediction/predictions_majority-vote_combined.txt"), emit: predictions_combined
    tuple val(meta), path("traitar_out/phenotype_prediction/predictions_majority-vote_phypat.txt"), emit: predictions_phypat
    tuple val(meta), path("traitar_out/phenotype_prediction/heatmap_*.png"), optional: true, emit: heatmaps
    tuple val(meta), path("traitar_out/pfam_annotation/*"), optional: true, emit: pfam_annotation
    tuple val(meta), path("traitar_out/gene_prediction/*"), optional: true, emit: gene_prediction
    tuple val("${task.process}"), val('traitar'), path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_type = task.ext.input_type ?: 'from_genes'
    // Validate input_type against allowed values
    if (!['from_genes', 'from_proteins', 'from_nucleotides'].contains(input_type)) {
        error("Invalid input_type: ${input_type}. Must be one of: 'from_genes', 'from_proteins', 'from_nucleotides'")
    }

    """
    # Create input directory with protein sequences
    mkdir -p input_dir
    cp ${proteins} input_dir/

    # Create sample file (tab-separated: filename, sample_name)
    echo -e "sample_file_name\\tsample_name" > samples.txt
    echo -e "${proteins.name}\\t${meta.id}" >> samples.txt

    # Run traitar phenotype analysis
    traitar phenotype \\
        input_dir \\
        samples.txt \\
        '${input_type}' \\
        traitar_out \\
        -c ${task.cpus} \\
        ${args}

    # Generate versions file
    cat > versions.yml <<-EOF
    "${task.process}":
        traitar: \$(traitar --version 2>&1 | grep -oE 'version.*' || echo 'unknown')
    EOF
    """

    stub:
    """
    mkdir -p traitar_out/phenotype_prediction
    mkdir -p traitar_out/pfam_annotation
    mkdir -p traitar_out/gene_prediction

    # Create phenotype prediction files with header and sample data
    cat > traitar_out/phenotype_prediction/predictions_majority-vote_combined.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    cat > traitar_out/phenotype_prediction/predictions_majority-vote_phypat.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    # Create a minimal PNG file (valid PNG header)
    printf '\\x89PNG\\r\\n\\x1a\\n' > traitar_out/phenotype_prediction/heatmap_combined.png

    # Create PFAM annotation with headers
    cat > traitar_out/pfam_annotation/${meta.id}.tsv <<-'EOF'
    protein_id\tpfam_id\tdomain_name\tevalue
    protein_001\tPF00001\tSample_domain_1\t1e-50
    protein_002\tPF00002\tSample_domain_2\t1e-45
    EOF

    # Create gene prediction with headers
    cat > traitar_out/gene_prediction/${meta.id}.tsv <<-'EOF'
    gene_id\tstart\tend\tstrand\tproduct
    gene_001\t1000\t2000\t+\tProtein_A
    gene_002\t2500\t3500\t-\tProtein_B
    EOF

    cat > versions.yml <<-EOF
    "${task.process}":
        traitar: 3.0.1
    EOF
    """
}
