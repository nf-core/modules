process LLAMACPPPYTHON_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c1/c16400973e9552f536824482b8bd1e02b48552cda1eee9921022e96b9ab8881c/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c1/c16400973e9552f536824482b8bd1e02b48552cda1eee9921022e96b9ab8881c/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/llama-cpp-python_llama.cpp:418f83c90cca5e9a' : 'community.wave.seqera.io/library/llama-cpp-python:0.3.28--0e4a0cb3e139c297') }"

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

    """
    echo ${args}
    """

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
