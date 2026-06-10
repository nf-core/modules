process LLAMACPPPYTHON_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'oras://community.wave.seqera.io/library/llama-cpp-python_llama.cpp:5fafadf87c4d2d13' : 'oras://community.wave.seqera.io/library/llama-cpp-python:0.3.28--aff015fe7d46e9ec') :
        (task.accelerator ? 'community.wave.seqera.io/library/llama-cpp-python_llama.cpp:0842a1eebc13be08' : 'community.wave.seqera.io/library/llama-cpp-python:0.3.28--c4e8129352f64626') }"

    input:
    tuple val(meta), path(prompt_file), path(gguf_model)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: output
    path "versions.yml", emit: versions_llama_cpp_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    template('llama-cpp-python.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        llama-cpp-python: \$(python3 -c 'import llama_cpp; print(llama_cpp.__version__)')
        cuda: no CUDA available
    END_VERSIONS
    """
}
