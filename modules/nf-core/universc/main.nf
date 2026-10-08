process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://wave.seqera.io/wt/31a6e5f54e99/wave/build:b9dbc80cb5dd84ee' :
        'wave.seqera.io/wt/31a6e5f54e99/wave/build:b9dbc80cb5dd84ee' }"

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
    # ------------------------------------------------------------
    # Detect Cell Ranger
    # ------------------------------------------------------------

    cr_version=\$(cellranger --version | sed 's/cellranger //')
    img_cr="/opt/cellranger-\${cr_version}"

    if [[ ! -d "\$img_cr" ]]; then
        echo "ERROR: Cell Ranger installation not found: \$img_cr" >&2
        exit 1
    fi


    # ------------------------------------------------------------
    # Create task-local Cell Ranger installation
    #
    # UniverSC modifies:
    #   lib/python
    #   mro/rna
    #   lib/python/cellranger/barcodes
    #
    # Everything else can remain in the immutable image.
    # ------------------------------------------------------------

    local_cr="\$PWD/.local-cellranger/cellranger-\${cr_version}"

    mkdir -p "\$local_cr"

    # Writable Cell Ranger Python code
    mkdir -p "\$local_cr/lib"
    cp -a "\${img_cr}/lib/python" "\$local_cr/lib/"

    # Writable Cell Ranger MRO
    mkdir -p "\$local_cr/mro"
    cp -a "\${img_cr}/mro/rna" "\$local_cr/mro/"

    # ------------------------------------------------------------
    # Symlink everything else from the image
    # ------------------------------------------------------------

    for item in "\${img_cr}"/*; do
        name=\$(basename "\$item")

        case "\$name" in
            cellranger|lib|mro|cellranger-cs)
                ;;
            *)
                ln -s "\$item" "\$local_cr/\$name"
                ;;
        esac
    done

    # Remaining lib contents
    for item in "\${img_cr}/lib"/*; do
        name=\$(basename "\$item")

        if [[ "\$name" != "python" ]]; then
            ln -s "\$item" "\$local_cr/lib/\$name"
        fi
    done

    # Remaining MRO contents
    for item in "\${img_cr}/mro"/*; do
        name=\$(basename "\$item")

        if [[ "\$name" != "rna" ]]; then
            ln -s "\$item" "\$local_cr/mro/\$name"
        fi
    done

    # Cell Ranger executable
    ln -s "\${img_cr}/cellranger" "\$local_cr/cellranger"

    # ------------------------------------------------------------
    # Create writable UniverSC installation
    # ------------------------------------------------------------

    mkdir -p "\$PWD/.local-universc"
    cp -a /opt/universc "\$PWD/.local-universc/"

    local_universc="\$PWD/.local-universc/universc"

    ln -s "\$local_universc/launch_universc.sh" "\$local_universc/universc"

    export PATH="\$local_cr:\$local_universc:\$PATH"

    # ------------------------------------------------------------
    # Run UniverSC
    # ------------------------------------------------------------

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
    # ------------------------------------------------------------
    # Detect Cell Ranger
    # ------------------------------------------------------------

    cr_version=\$(cellranger --version | sed 's/cellranger //')
    img_cr="/opt/cellranger-\${cr_version}"

    if [[ ! -d "\$img_cr" ]]; then
        echo "ERROR: Cell Ranger installation not found: \$img_cr" >&2
        exit 1
    fi


    # ------------------------------------------------------------
    # Create task-local Cell Ranger installation
    #
    # UniverSC modifies:
    #   lib/python
    #   mro/rna
    #   lib/python/cellranger/barcodes
    #
    # Everything else can remain in the immutable image.
    # ------------------------------------------------------------

    local_cr="\$PWD/.local-cellranger/cellranger-\${cr_version}"

    mkdir -p "\$local_cr"

    # Writable Cell Ranger Python code
    mkdir -p "\$local_cr/lib"
    cp -a "\${img_cr}/lib/python" "\$local_cr/lib/"

    # Writable Cell Ranger MRO
    mkdir -p "\$local_cr/mro"
    cp -a "\${img_cr}/mro/rna" "\$local_cr/mro/"

    # ------------------------------------------------------------
    # Symlink everything else from the image
    # ------------------------------------------------------------

    for item in "\${img_cr}"/*; do
        name=\$(basename "\$item")

        case "\$name" in
            cellranger|lib|mro|cellranger-cs)
                ;;
            *)
                ln -s "\$item" "\$local_cr/\$name"
                ;;
        esac
    done

    # Remaining lib contents
    for item in "\${img_cr}/lib"/*; do
        name=\$(basename "\$item")

        if [[ "\$name" != "python" ]]; then
            ln -s "\$item" "\$local_cr/lib/\$name"
        fi
    done

    # Remaining MRO contents
    for item in "\${img_cr}/mro"/*; do
        name=\$(basename "\$item")

        if [[ "\$name" != "rna" ]]; then
            ln -s "\$item" "\$local_cr/mro/\$name"
        fi
    done

    # Cell Ranger executable
    ln -s "\${img_cr}/cellranger" "\$local_cr/cellranger"

    # ------------------------------------------------------------
    # Create writable UniverSC installation
    # ------------------------------------------------------------

    mkdir -p "\$PWD/.local-universc"
    cp -a /opt/universc "\$PWD/.local-universc/"

    local_universc="\$PWD/.local-universc/universc"

    ln -s "\$local_universc/launch_universc.sh" "\$local_universc/universc"

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
