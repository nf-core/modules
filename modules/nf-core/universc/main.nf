process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/nf-core/universc:1.2.5.1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val(technology)

    output:
    tuple val(meta), path("${prefix}/outs/*"), emit: outs
    tuple val("${task.process}"), val('cellranger'), eval("cellranger 2>&1 | sed '/^cellranger/!d;s/cellranger  (//;s/)//'"), emit: versions_cellranger, topic: versions
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
    cr_version=\$(cellranger 2>&1 | sed '/^cellranger/!d;s/cellranger  (//;s/)//')
    img_cr="/cellranger-\$cr_version"
    img_cs="\$img_cr/cellranger-cs/\$cr_version"

    local_cr="\$PWD/.local-cellranger/cellranger-\$cr_version"
    local_cs="\$local_cr/cellranger-cs/\$cr_version"

    # Create local Cell Ranger directory structure
    mkdir -p "\$local_cs/lib/python" "\$local_cs/mro"

    # Top-level Cell Ranger files/directories
    ln -s \$img_cr/cellranger-tiny-fastq "\$local_cr/cellranger-tiny-fastq"
    ln -s \$img_cr/cellranger-tiny-ref  "\$local_cr/cellranger-tiny-ref"

    # Cell Ranger CS: everything except lib and mro is symlinked
    for item in \$img_cs/*; do
        name=\$(basename "\$item")
        if [[ "\$name" != "lib" && "\$name" != "mro"  ]]; then
            ln -s "\$item" "\$local_cs/\$name"
        fi
    done

    # lib: everything except python is symlinked
    for item in \$img_cs/lib/*; do
        name=\$(basename "\$item")
        if [[ "\$name" != "python" ]]; then
            ln -s "\$item" "\$local_cs/lib/\$name"
        fi
    done

    # Copy python and mro folders (~191 MB, mostly barcodes) as modified by universc
    cp -a \$img_cs/lib/python "\$local_cs/lib"
    cp -a \$img_cs/mro \$local_cs

    # Symlink cellranger bin
    ln -s "\$local_cs/bin/cellranger" "\$local_cr/cellranger"

    # UNIVERSC needs its installation directory to be writable
    mkdir -p "\$PWD/.local-universc"
    cp -a /universc "\$PWD/.local-universc/"

    local_universc="\$PWD/.local-universc/universc"

    # Fix UNIVERSC launcher symlink
    rm "\$local_universc/universc"
    ln -s "\$local_universc/launch_universc.sh" "\$local_universc/universc"

    # Fix malformed [[ syntax in UNIVERSC
    sed -i 's/"\$technology" == "vasa-drop"/ "\$technology" == "vasa-drop"/' "\$local_universc/universc"

    export PATH="\$local_cr:\$local_universc:\$PATH"

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
    cr_version=\$(cellranger 2>&1 | sed '/^cellranger/!d;s/cellranger  (//;s/)//')
    img_cr="/cellranger-\$cr_version"
    img_cs="\$img_cr/cellranger-cs/\$cr_version"

    local_cr="\$PWD/.local-cellranger/cellranger-\$cr_version"
    local_cs="\$local_cr/cellranger-cs/\$cr_version"

    # Create local Cell Ranger directory structure
    mkdir -p "\$local_cs/lib/python" "\$local_cs/mro"

    # Top-level Cell Ranger files/directories
    ln -s \$img_cr/cellranger-tiny-fastq "\$local_cr/cellranger-tiny-fastq"
    ln -s \$img_cr/cellranger-tiny-ref  "\$local_cr/cellranger-tiny-ref"

    # Cell Ranger CS: everything except lib and mro is symlinked
    for item in \$img_cs/*; do
        name=\$(basename "\$item")
        if [[ "\$name" != "lib" && "\$name" != "mro"  ]]; then
            ln -s "\$item" "\$local_cs/\$name"
        fi
    done

    # lib: everything except python is symlinked
    for item in \$img_cs/lib/*; do
        name=\$(basename "\$item")
        if [[ "\$name" != "python" ]]; then
            ln -s "\$item" "\$local_cs/lib/\$name"
        fi
    done

    # Copy python and mro folders (~191 MB, mostly barcodes) as modified by universc
    cp -a \$img_cs/lib/python "\$local_cs/lib"
    cp -a \$img_cs/mro \$local_cs

    # Symlink cellranger bin
    ln -s "\$local_cs/bin/cellranger" "\$local_cr/cellranger"

    # UNIVERSC needs its installation directory to be writable
    mkdir -p "\$PWD/.local-universc"
    cp -a /universc "\$PWD/.local-universc/"

    local_universc="\$PWD/.local-universc/universc"

    # Fix UNIVERSC launcher symlink
    rm "\$local_universc/universc"
    ln -s "\$local_universc/launch_universc.sh" "\$local_universc/universc"

    # Fix malformed [[ syntax in UNIVERSC
    sed -i 's/"\$technology" == "vasa-drop"/ "\$technology" == "vasa-drop"/' "\$local_universc/universc"

    export PATH="\$local_cr:\$local_universc:\$PATH"

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
