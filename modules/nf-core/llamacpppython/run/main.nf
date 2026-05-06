process LLAMACPPPYTHON_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f5c3807fc8d9a8c28cef653f377cec46d26097ad157571001316dd8ca196fb4/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3c/3c1e6f266b077ae89b341bf921b59a3d5d382104c3775725a814bd68f064541e/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_pip_cuda-version_cuda-runtime:202c39541466fecd' : 'community.wave.seqera.io/library/llama-cpp-python:0.3.16--b351398cd0ea7fc5') }"

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
