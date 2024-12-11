
process IGV_JS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(alignment), path(index)

    output:
    tuple val(meta), path("*_genome-browser.html") , emit: browser
    tuple val(meta), path(alignment)               , emit: align_files
    tuple val(meta), path(index)                   , emit: index_files
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def browser_args = task.ext.args ?: ''
    def track_args = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat  <<IGV > ${prefix}_genome-browser.html
    <html>
        <head>
            <script src="https://cdn.jsdelivr.net/npm/igv@2.13.9/dist/igv.min.js"></script> <!-- https://github.com/igvteam/igv.js/ -->
        </head>
        <body>

        <div id="igv-div"></div>
        <script>
            var igvDiv = document.getElementById("igv-div");
            var options =
            {
                $browser_args,
                tracks: [
                    {
                        "name": "$prefix",
                        "url": "$alignment",
                        "indexURL": "$index",
                        $track_args
                    }
                ]
            };
            igv.removeAllBrowsers();
            igv.createBrowser(igvDiv, options)
                .then(function(browser){
                    igv.browser = browser;
                })
        </script>

        </body>
    </html>
    IGV

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genome-browser.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
