include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process METABAT2_METABAT2 {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::metabat2=2.15" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/metabat2:2.15--h986a166_1"
    } else {
        container "quay.io/biocontainers/metabat2:2.15--h986a166_1"
    }

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    tuple val(meta), path("metabat2/*.fa")   , optional:true , emit: fasta
    tuple val(meta), path("metabat2/$prefix"), optional:true , emit: membership
    path "versions.yml"                                      , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def depth_file = depth.size() == 0 ? "" : "-a $depth"
    """
    metabat2 \\
        $options.args \\
        -i $fasta \\
        $depth_file \\
        -t $task.cpus \\
        --saveCls \\
        -o metabat2/${prefix}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """
}
