process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://wave.seqera.io/wt/9798dc70dd2e/wave/build:ac5e4fbf01375a84' :
        'wave.seqera.io/wt/9798dc70dd2e/wave/build:ac5e4fbf01375a84' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val(technology)

    output:
    tuple val(meta), path("${prefix}/outs/*"), emit: outs
    tuple val("${task.process}"), val('cellranger'), eval("cellranger --version | sed 's/cellranger //'"), emit: versions_cellranger, topic: versions
    tuple val("${task.process}"), val('universc'), eval("universc --version | sed -n 's/launch_universc.sh version //p'"), emit: versions_universc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args        = task.ext.args   ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--file $reads" : "-R1 ${reads[0]} -R2 ${reads[1]}"

    def reference_name = reference.name
    """
    export PATH="/opt/cellranger-10.1.0:\$PATH"
    sed -i 's/"\$technology" == "vasa-drop"/ "\$technology" == "vasa-drop"/' "/opt/universc/launch_universc.sh"
    sed -i 's/\$index1(\\[0\\])/"\${index1[0]}"/g'  "/opt/universc/launch_universc.sh"
    sed -i 's/\${index2}(\\[0\\])/"\${index2[0]}"/g'  "/opt/universc/launch_universc.sh"
    sed -i 's/\$index2(\\[0\\])/"\${index2[0]}"/g'  "/opt/universc/launch_universc.sh"
    sed -i 's/bam="--no-bam"/bam="--create-bam=false"/' "/opt/universc/launch_universc.sh"
    sed -i 's/bam=""/bam="--create-bam=true"/' "/opt/universc/launch_universc.sh"

    sed -i '2640c\
        if false; then
    ' /opt/universc/launch_universc.sh

    sed -i '2690c\
        elif false; then
    ' /opt/universc/launch_universc.sh
    sed -n '2637,2643p' /opt/universc/launch_universc.sh

    export PYTHON_EGG_CACHE=\$(pwd)/.cache
    universc \\
        --id ${prefix} \\
        ${input_reads} \\
        --technology ${technology} \\
        --reference ${reference_name} \\
        --jobmode "local" \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        --per-cell-data \\
        ${args}

    # save log files
    echo !! > ${prefix}/outs/_invocation
    """


    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export PATH="/opt/cellranger-10.1.0:\$PATH"
    sed -i 's/"\$technology" == "vasa-drop"/ "\$technology" == "vasa-drop"/' "/opt/universc/launch_universc.sh"
    sed -i 's/\$index1(\\[0\\])/"\${index1[0]}"/g'  "/opt/universc/launch_universc.sh"
    sed -i 's/\${index2}(\\[0\\])/"\${index2[0]}"/g'  "/opt/universc/launch_universc.sh"
    sed -i 's/\$index2(\\[0\\])/"\${index2[0]}"/g'  "/opt/universc/launch_universc.sh"
    sed -i 's/bam="--no-bam"/bam="--create-bam=false"/' "/opt/universc/launch_universc.sh"
    sed -i 's/bam=""/bam="--create-bam=true"/' "/opt/universc/launch_universc.sh"

    sed -i '2640c\
        if false; then
    ' /opt/universc/launch_universc.sh

    sed -i '2690c\
        elif false; then
    ' /opt/universc/launch_universc.sh
    sed -n '2637,2643p' /opt/universc/launch_universc.sh

    mkdir -p ${prefix}/outs/
    cd ${prefix}/outs/

    touch _invocation

    touch basic_stats.txt
    touch metrics_summary.csv
    touch molecule_info.h5
    touch possorted_genome_bam.bam
    touch possorted_genome_bam.bam.bai
    touch web_summary.html

    mkdir -p filtered_feature_bc_matrix
    touch filtered_feature_bc_matrix.h5
    echo "" | gzip > filtered_feature_bc_matrix/barcodes.tsv.gz
    echo "" | gzip > filtered_feature_bc_matrix/features.tsv.gz
    echo "" | gzip > filtered_feature_bc_matrix/matrix.mtx.gz

    mkdir -p raw_feature_bc_matrix
    touch raw_feature_bc_matrix.h5
    echo "" | gzip > raw_feature_bc_matrix/barcodes.tsv.gz
    echo "" | gzip > raw_feature_bc_matrix/features.tsv.gz
    echo "" | gzip > raw_feature_bc_matrix/matrix.mtx.gz

    cd ../..
    """
}
