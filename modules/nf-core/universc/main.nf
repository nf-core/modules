process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'

    container "nf-core/universc:1.2.5.1"
    containerOptions {
        ['singularity', 'apptainer'].contains(workflow.containerEngine)
            ? "-B /var/tmp --writable-tmpfs"
            : workflow.containerEngine == 'docker'
                ? "--privileged"
                : workflow.containerEngine == 'podman'
                    ? "--runtime crun --userns=keep-id --systemd=always"
                    : ''
    }

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("sample-${meta.id}/outs/*"), emit: outs
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args        = task.ext.args   ?: ''
    def input_reads = meta.single_end ? "--file $reads" : "-R1 ${reads[0]} -R2 ${reads[1]}"

    def sample_arg     = meta.samples.unique().join(",")
    def reference_name = reference.name
    """
    export PYTHON_EGG_CACHE=\$(pwd)/.cache
    universc \\
        --id 'sample-${meta.id}' \\
        ${input_reads} \\
        --technology '${meta.technology}' \\
        --chemistry '${meta.chemistry}' \\
        --reference ${reference_name} \\
        --description ${sample_arg} \\
        --jobmode "local" \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        --per-cell-data \\
        ${args} 1> _log 2> _err

    # save log files
    echo !! > sample-${meta.id}/outs/_invocation
    cp _log sample-${meta.id}/outs/_log
    cp _err sample-${meta.id}/outs/_err

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger:  \$(echo \$(cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g' ))
        universc:  \$(echo \$(bash /universc/launch_universc.sh --version | grep version | grep universc  | sed 's/^.* //g' ))
    END_VERSIONS
    """


    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    mkdir -p sample-${meta.id}/outs/
    cd sample-${meta.id}/outs/

    touch _invocation
    touch _log
    touch _err

    touch basic_stats.txt
    touch metrics_summary.csv
    touch molecule_info.h5
    touch possorted_genome_bam.bam
    touch possorted_genome_bam.bam.bai
    touch web_summary.html

    mkdir -p filtered_feature_bc_matrix
    touch filtered_feature_bc_matrix.h5
    echo | gzip > filtered_feature_bc_matrix/barcodes.tsv.gz
    echo | gzip > filtered_feature_bc_matrix/features.tsv.gz
    echo | gzip > filtered_feature_bc_matrix/matrix.mtx.gz

    mkdir -p raw_feature_bc_matrix
    touch raw_feature_bc_matrix.h5
    echo | gzip > raw_feature_bc_matrix/barcodes.tsv.gz
    echo | gzip > raw_feature_bc_matrix/features.tsv.gz
    echo | gzip > raw_feature_bc_matrix/matrix.mtx.gz

    cd ../..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger:  \$(echo \$(cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g' ))
        universc:  \$(echo \$(bash /universc/launch_universc.sh --version | grep version | grep universc | sed 's/^.* //g' ))
    END_VERSIONS
    """
}
