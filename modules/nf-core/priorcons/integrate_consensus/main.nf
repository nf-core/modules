
process PRIORCONS_INTEGRATE_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/priorcons:0.1.0--pyhdfd78af_0' :
    'biocontainers/priorcons:0.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_aln), val(ref_id), path(prior_parquet)

    output:
    tuple val(meta), path("${meta.id}_output"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    priorcons integrate-consensus \\
        --input $input_aln \\
        --ref $ref_id \\
        --prior $prior_parquet \\
        --output_dir ${prefix}_output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        priorcons: \$(priorcons --version 2>&1 | head -n1 || echo "v0.1.0")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    mkdir -p ${prefix}_output
    
    # Create stub out files
    echo ">${prefix}_integrated_consensus" > ${prefix}_output/${prefix}-INTEGRATED.fasta
    echo "ATCGATCGATCG" >> ${prefix}_output/${prefix}-INTEGRATED.fasta
    
    # Create windows_trace.csv stub
    cat <<-EOF > ${prefix}_output/windows_trace.csv
    window_start,window_end,score,selected_sequence
    1,100,0.95,mapping_consensus
    101,200,0.87,abacas_consensus
    EOF
    
    # Create qc.json stub
    cat <<-EOF > ${prefix}_output/qc.json
    {
        "total_windows": 150,
        "mapping_selected": 95,
        "abacas_selected": 55,
        "average_score": 0.91
    }
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        priorcons: "v0.1.0"
    END_VERSIONS
    """
}