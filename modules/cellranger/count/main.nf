// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CELLRANGER_COUNT {
    tag "$meta.gem"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "nfcore/cellranger:6.0.2"

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    path "versions.yml"               , emit: versions
    path("sample-${meta.gem}/outs/*") , emit: outs

    script:
    def sample_arg = meta.samples.unique().join(",")
    def reference_name = reference.name
    """
    cellranger count --id='sample-${meta.gem}' \
        --fastqs=. \
        --transcriptome=${reference_name} \
        --sample=${sample_arg} \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()} \
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "sample-${meta.gem}/outs/"
    touch sample-${meta.gem}/outs/fake_file.txt
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
