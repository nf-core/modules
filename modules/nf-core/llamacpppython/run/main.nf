process LLAMACPPPYTHON_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c3240996011a266ad1bf60a9591b995f72e17d80d28827f319ae8f4059c720c3/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/185e274b6759747466cb5b0a1f7c0c7a462e790c2b783b225d260bbd8526644c/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/llama-cpp-python_llama.cpp_cuda-version:8eda7e45fafc0673' : 'community.wave.seqera.io/library/llama-cpp-python_llama.cpp:211706178075740f') }"

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
