process SIGPROFILER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b6141a7b8f0674ac604d90eb1306a731da24a734:daf9213409b038023df2a741f058d1bfd66b3c4c-0':
        'biocontainers/mulled-v2-b6141a7b8f0674ac604d90eb1306a731da24a734:daf9213409b038023df2a741f058d1bfd66b3c4c-0' }"

    input:
    tuple val(meta), path(tsv_list, stageAs: '*.tsv')
    
    output:
    tuple val(meta), path("results/*"),    emit: sigprofiler_results
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.py"
   

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p results
    echo "stub output" > results/dummy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigprofilermatrixGenerator: \$(python3 -c "import SigProfilerMatrixGenerator as matGen; print(matGen.__version__)")
        sigprofilerextractor: \$(python3 -c "import SigProfilerExtractor as sig; print(sig.__version__)")
    END_VERSIONS
    """
}
