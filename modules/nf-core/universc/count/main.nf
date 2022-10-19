process UNIVERSC_COUNT {
    tag "$meta.id"
    label 'process_medium'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
        conda (params.enable_conda ? "hcc::cellranger=3.0.2" : null)
    }
    container "tomkellygenetics/universc:1.2.4"
    containerOptions = "--user root"

    input:
    tuple val(meta), path(reads)
    path  reference

    output:
    tuple val(meta), path("sample-${meta.id}/outs/*"), emit: outs
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_arg = meta.samples.unique().join(",")
    def reference_name = reference.name
    """
    cellranger \\
        count  \\
        --id='sample-${meta.id}' \\
        --fastqs=. \\
        --transcriptome=${reference_name} \\
        --sample=${sample_arg} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger:  \$(echo \$(cellranger count --version 2>&1 | head -n 2 | tail -n 1 | sed 's/^.* //g' | sed 's/(//g' | sed 's/)//g' ))
    END_VERSIONS
    """
}
