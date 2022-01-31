//def VERSION=0.0.1

process CUSTOM_METHODSDESCRIPTION {
    label 'process_low'

    output:
    path "*_mqc.html"   , emit: mqc_html
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def doi = workflow.manifest.doi ? " (${workflow.manifest.doi})" : ""

    """
    cat <<-METHOD_DESCRIPTION > method_description_mqc.html
    <!--
    id: 'methods-description'
    section_name: 'Methods Description'
    -->
    <p>
    Data was processed using the nf-core (Ewels et al. 2020) pipeline ${workflow.manifest.name} v${workflow.manifest.version}$doi.

    The pipeline was executed with Nextflow (v${nextflow.version}; Di Tommaso et al. 2017) with the following command:</p>

    <pre><code>
    ${workflow.commandLine}
    </code></pre>

    <b>References</b>
    <ul>
        <li>Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. https://doi.org/10.1038/nbt.3820 </li>
        <li>Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38(3), 276–278. https://doi.org/10.1038/s41587-020-0439-x </li>
    </ul>

    <blockquote>
    <b>Notes:</b><br>
    <ul>
        <li>If present, make sure to replace the pipeline DOI with correct reference information.</li>
        <li>The command above does not include parameters contained in any configs that may have been used. Ensure the config file is also uploaded with your publication.</li>
        <li>You should also cite all software used within this run. Check the 'Software Versions' of this report to get version information.</li>
    </ul>
    </blockquote>
    METHOD_DESCRIPTION

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        custom: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
