// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process TRAITAR {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::traitar=3.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://a_gihawi/traitar3/traitar3' :
        'quay.io/biocontainers/traitar:3.0.1--pyh7cba7a3_0' }"

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
        ${input_type} \\
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
