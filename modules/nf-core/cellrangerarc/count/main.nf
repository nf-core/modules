process CELLRANGERARC_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "heylf/cellranger-arc:2.0.2"

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use docker or singularity containers."
    }

    input:
    tuple val(meta), path(reads)
    path  lib_csv // A script to generate the lib.csv automatically you can find at nf-core/scrnaseq
    path  reference

    output:
    tuple val(meta), path("${meta.id}/outs/*"), emit: outs
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def lib_csv_name = lib_csv.name
    def reference_name = reference.name
    """
    # The following ugly three commands (mkdir, mv, awk) are required because cellranger-arc only deals with abolsute paths
    if [ ! -d "fastqs" ]; then
        mkdir fastqs
    fi

    mv *.fastq.gz fastqs/

    # This adds the tmp/fastqs to the lib.csv
    cat $lib_csv_name | awk '{if (\$0 !~ /^fastqs/) {print "'"\$(readlink -f fastqs)"'"\$0} else {print \$0}}' > tmp.csv && mv tmp.csv $lib_csv_name

    cellranger-arc \\
        count \\
        --id='${meta.id}' \\
        --libraries=$lib_csv_name \\
        --reference=$reference_name \\
        --localcores=$task.cpus \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangerarc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}/outs/"
    touch ${meta.id}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangerarc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
