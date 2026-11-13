nextflow.enable.types = true

process CUSTOM_WRITEFILE {
    tag name

    input:
    record(
        name: String,
        items: List<String>,
        newLine: Boolean?
    )

    output:
    file(name)

    exec:
    def path = task.workDir.resolve(name)
    path.delete()
    items.each { item ->
        path << item
        if (newLine) {
            path << '\n'
        }
    }
}
