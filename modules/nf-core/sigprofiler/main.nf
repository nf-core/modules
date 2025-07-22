process SIGPROFILER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b6141a7b8f0674ac604d90eb1306a731da24a734:daf9213409b038023df2a741f058d1bfd66b3c4c-0':
        'biocontainers/mulled-v2-b6141a7b8f0674ac604d90eb1306a731da24a734:daf9213409b038023df2a741f058d1bfd66b3c4c-0' }"

    input:
    tuple val(meta), path(tsv_list, stageAs: '*.tsv')
    path(genome_installed_path)  //optional
    
    output:
    tuple val(meta), path("results/SBS96/*")  , emit: sbs_results
    tuple val(meta), path("results/ID83/*")   , emit: id_results, optional: true
    tuple val(meta), path("results/DBS78/*")  , emit: dbs_results, optional: true 
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.py"

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p results/SBS96/SBS96/
    touch results/SBS96/SBS96/Samples.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
        sigprofilermatrixGenerator: \$(python3 -c "import SigProfilerMatrixGenerator as matGen; print(matGen.__version__)")
        sigprofilerextractor: \$(python3 -c "import SigProfilerExtractor as sig; print(sig.__version__)")
    END_VERSIONS
    """
}
