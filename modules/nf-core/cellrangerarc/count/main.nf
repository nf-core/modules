process CELLRANGERARC_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger-arc:2.0.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), val(sample_type), val(sub_sample), path(reads, stageAs: "fastqs/*")
    path  reference

    output:
    tuple val(meta), path("${meta.id}/outs/**"), emit: outs
    path("${meta.id}_lib.csv")                 , emit: lib
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_name = reference.name
    def sample_types = sample_type.join(",")
    def sample_names = sub_sample.join(",")
    def lib_csv = meta.id + "_lib.csv"

    """
    fastq_folder=\$(readlink -f fastqs)

    python3 <<CODE

    sample_types = "${sample_types}".split(",")
    sample_names = "${sample_names}".split(",")
    unique_samples_names = set(sample_names)

    lib_csv = open("${lib_csv}", "w")
    lib_csv.write("fastqs,sample,library_type")

    for i in range(0, len(sample_types)):
        if sample_names[i] in unique_samples_names:
            unique_samples_names.remove(
                sample_names[i]
            )  # this has to be done to account for different Lane files (e.g., L002)
            if sample_types[i] == "gex":
                lib_csv.write("\\n{},{},{}".format("\${fastq_folder}", sample_names[i], "Gene Expression"))
            else:
                lib_csv.write("\\n{},{},{}".format("\${fastq_folder}", sample_names[i], "Chromatin Accessibility"))

    lib_csv.close()

    print("Wrote lib.csv file to {}".format("${lib_csv}"))

    CODE

    cellranger-arc \\
        count \\
        --id='${meta.id}' \\
        --libraries=$lib_csv \\
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
