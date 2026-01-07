process TRAITAR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/hmmer_prodigal_pandas_pip_pruned:a83f0296374a52e6"

    input:
    tuple val(meta), path(proteins)

    output:
    tuple val(meta), path("traitar_out/phenotype_prediction/predictions_majority-vote_combined.txt"), emit: predictions_combined
    tuple val(meta), path("traitar_out/phenotype_prediction/predictions_single-votes_combined.txt"), emit: predictions_single_votes
    tuple val(meta), path("traitar_out/phenotype_prediction/predictions_flat_*.txt"), optional: true, emit: predictions_flat
    tuple val(meta), path("traitar_out/predictions_*.txt"), optional: true, emit: predictions_raw
    tuple val(meta), path("traitar_out/annotation/pfam/"), optional: true, emit: pfam_annotation
    tuple val(meta), path("traitar_out/gene_prediction/"), optional: true, emit: gene_prediction
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

    # Handle gzipped files by decompressing
    if [[ "${proteins}" == *.gz ]]; then
        gunzip -c "${proteins}" > "input_dir/\$(basename '${proteins}' .gz)"
        INPUT_FILE="\$(basename '${proteins}' .gz)"
    else
        cp "${proteins}" input_dir/
        INPUT_FILE="\$(basename '${proteins}')"
    fi

    # Create sample file (tab-separated: filename, sample_name)
    echo -e "sample_file_name\\tsample_name" > samples.txt
    echo -e "\${INPUT_FILE}\\t${meta.id}" >> samples.txt

    # Download PFAM data for traitar annotation
    mkdir -p pfam_data
    traitar pfam pfam_data

    # Run traitar phenotype analysis
    traitar phenotype \\
        pfam_data \\
        input_dir \\
        samples.txt \\
        '${input_type}' \\
        traitar_out \\
        -c ${task.cpus} \\
        --overwrite \\
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
    mkdir -p traitar_out/annotation/pfam
    mkdir -p traitar_out/gene_prediction

    # Create phenotype prediction files with header and sample data
    cat > traitar_out/phenotype_prediction/predictions_majority-vote_combined.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    cat > traitar_out/phenotype_prediction/predictions_single-votes_combined.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    # Create PFAM annotation with headers
    cat > traitar_out/annotation/pfam/${meta.id}_filtered_best.dat <<-'EOF'
    protein_id\tpfam_id\tdomain_name\tevalue
    protein_001\tPF00001\tSample_domain_1\t1e-50
    protein_002\tPF00002\tSample_domain_2\t1e-45
    EOF

    cat > traitar_out/annotation/pfam/${meta.id}_domtblout.dat <<-'EOF'
    target name\taccession\tquery name\taccession\tE-value\tscore
    PF00001.10\tPF00001.10\tprotein_001\t-\t1e-50\t100.0
    PF00002.10\tPF00002.10\tprotein_002\t-\t1e-45\t95.0
    EOF

    cat > traitar_out/annotation/pfam/summary.dat <<-'EOF'
    protein_id\ttotal_domains\tstatus
    protein_001\t1\tannotated
    protein_002\t1\tannotated
    EOF

    # Create raw prediction files
    cat > traitar_out/predictions_raw.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0.5\t0.8\t0.2
    EOF

    cat > traitar_out/predictions_majority-vote.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    cat > traitar_out/predictions_single-votes.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    cat > traitar_out/predictions_conservative-vote.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    # Create flat format prediction files
    cat > traitar_out/phenotype_prediction/predictions_flat_majority-votes_combined.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF

    cat > traitar_out/phenotype_prediction/predictions_flat_single-votes_combined.txt <<-'EOF'
    sample_name\tphenotype_1\tphenotype_2\tphenotype_3
    ${meta.id}\t0\t1\t0
    EOF
    cat > traitar_out/gene_prediction/${meta.id}.faa <<-'EOF'
    >gene_001
    MIVLTGGNAAGLLAALLGAPVVAIAGKGTDKKSAAA
    >gene_002
    MGSTNGAKPQSLKQALVQAPQGVSGKEVVGSPPEAT
    EOF

    cat > traitar_out/gene_prediction/${meta.id}.gff <<-'EOF'
    ##gff-version 3
    sequence\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes
    test_sample\tprodigal\tCDS\t1000\t2000\t.\t+\t0\tID=gene_001
    test_sample\tprodigal\tCDS\t2500\t3500\t.\t-\t0\tID=gene_002
    EOF

    cat > versions.yml <<-EOF
    "${task.process}":
        traitar: 3.0.1
    EOF
    """
}
