String dedent(String input) {
    def lines = input.tokenize('\n')
    def minIndent = Integer.MAX_VALUE

    for (line in lines) {
        if (line.trim()) {
            def indent = line.indexOf(line.trim())
            if (indent < minIndent) {
                minIndent = indent
            }
        }
    }

    def dedentedLines = []
    for (line in lines) {
        dedentedLines << line.substring(minIndent)
    }

    return dedentedLines.join('\n');
}


String _make_versions_html(def versions) {
    def html = [
        dedent(
            """\\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
            """
        )
    ]

    versions.each { process, tmp_versions ->
        html << "<tbody>"
        tmp_versions.entrySet().sort { a, b -> a.key <=> b.key }.eachWithIndex { entry, i ->
            def tool = entry.key
            def version = entry.value
            html << 
                dedent(
                    """\\
                    <tr>
                        <td><samp>${i == 0 ? process : ''}</samp></td>
                        <td><samp>${tool}</samp></td>
                        <td><samp>${version}</samp></td>
                    </tr>
                    """
                )
            
        html << "</tbody>"
        }
    }
    html << "</table>"
    return html.join('\n')
}

def main() {
    def versionsThisModule = [:]
    versionsThisModule[task.process] = [
        "groovy": System.getProperty("groovy.version"),
        "yaml": org.yaml.snakeyaml.Yaml.class.getDeclaredField('version').get(null)
    ]

    @Grab('org.yaml:snakeyaml:1.29') // Add SnakeYAML library dependency
    import org.yaml.snakeyaml.Yaml
    import org.yaml.snakeyaml.Loader
    def versionsFile = new File(versions) // Assuming "versions" is the file path
    def versionsByProcess
    try {
        def yaml = new Yaml(new Loader(LoaderOptions.builder().allowDuplicateKeys(false).build()))
        def versionsByModule = yaml.load(versionsFile.text)
        versionsByProcess = versionsByModule + versionsThisModule
    } catch (Exception e) {
        println("Error: ${e.message}")
    }

    versionByModule = [:]
    versionsByProcess.each { process, processVersions -> 
        def module = process.tokenize(":")[-1]
        try {
            if (versionsByModule[module] != processVersions) {
                throw new AssertionError("We assume that software versions are the same between all modules. If you see this error-message it means you discovered an edge-case and should open an issue in nf-core/tools.")
            }
        }
        catch (MissingPropertyException) {
            versionByModule[module] = processVersions
        }
    }

    versionByModule["Workflow"] = [
        "Nextflow": workflow.nextflow.version,
        workflow.manifest.name: workflow.manifest.version,
    ]

    def versionsMqc = [
        "id": "software_versions",
        "section_name": "${workflow.manifest.name} Software Versions",
        "section_href": "https://github.com/${workflow.manifest.name}",
        "plot_type": "html",
        "description": "are collected at run time from the software output.",
        "data": _make_versions_html(versionByModule),
    ]

    @Grab('org.yaml:snakeyaml:1.29') // Add SnakeYAML library dependency
    import org.yaml.snakeyaml.Yaml

    // Write to "software_versions.yml"
    def softwareVersionsYaml = new Yaml()
    def softwareVersionsFile = new File("software_versions.yml")
    softwareVersionsFile.text = softwareVersionsYaml.dump(versionsByModule)
    // Write to "software_versions_mqc.yml"
    def softwareVersionsMqcYaml = new Yaml()
    def softwareVersionsMqcFile = new File("software_versions_mqc.yml")
    softwareVersionsMqcFile.text = softwareVersionsMqcYaml.dump(versionsMqc)
    // Write to "versions.yml"
    def versionsYaml = new Yaml()
    def versionsFile = new File("versions.yml")
    versionsFile.text = versionsYaml.dump(versionsThisModule)

}

main()