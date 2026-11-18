nextflow.enable.types = true

process CUSTOM_WRITEFILE {
    tag "${id}"

    ext [
        prefix: { output_prefix ?: "${id}" }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        output_prefix: String?,
        extension: String,
        items: List<String>,
        newLine: Boolean?
    )

    output:
    record(
        id: id,
        meta: meta,
        file: file("${task.ext.prefix}.${extension}")
    )

    exec:
    task.workDir.resolve("${task.ext.prefix}.${extension}").text = newLine ? items.collect { item -> "${item}\n" }.join() : items.join()
}
