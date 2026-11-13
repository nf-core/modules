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
    def path = task.workDir.resolve("${task.ext.prefix}.${extension}")
    path.delete()
    items.each { item ->
        path << item
        if (newLine) {
            path << '\n'
        }
    }
}
