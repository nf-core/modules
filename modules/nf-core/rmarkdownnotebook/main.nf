process RMARKDOWNNOTEBOOK {
    tag "$meta.id"
    label 'process_low'

    //NB: You likely want to override this with a container containing all required
    //dependencies for your analysis. The container at least needs to contain the
    //yaml and rmarkdown R packages.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' :
        'quay.io/biocontainers/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' }"

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files

    output:
    tuple val(meta), path("*.html")             , emit: report
    tuple val(meta), path("*.parameterised.Rmd"), emit: parameterised_notebook, optional: true
    tuple val(meta), path("artifacts/*")        , emit: artifacts, optional: true
    tuple val(meta), path("session_info.log")   , emit: session_info
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    notebook_parameters = [
        meta: meta,
        cpus: task.cpus,
    ] + (parameters ?: [:])

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def yamlBuilder = new groovy.yaml.YamlBuilder()
    yamlBuilder.call(notebook_parameters)
    def yaml_content = yamlBuilder.toString().tokenize('\n').join("\n    ")
    """
    # Dump parameters to yaml file
    cat <<- END_YAML_PARAMS > params.yml
    ${yaml_content}
    END_YAML_PARAMS

    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="$task.cpus"
    export OPENBLAS_NUM_THREADS="$task.cpus"
    export OMP_NUM_THREADS="$task.cpus"

    # Work around  https://github.com/rstudio/rmarkdown/issues/1508
    # If the symbolic link is not replaced by a physical file
    # output- and temporary files will be written to the original directory.
    mv "${notebook}" "${notebook}.orig"
    cp -L "${notebook}.orig" "${prefix}.Rmd"

    # Render notebook
    Rscript - <<EOF
    params = yaml::read_yaml('params.yml')

    # Instead of rendering with params, produce a version of the R
    # markdown with param definitions set, so the notebook itself can
    # be reused
    rmd_content <- readLines('${prefix}.Rmd')

    # Extract YAML content between the first two '---'
    start_idx <- which(rmd_content == "---")[1]
    end_idx <- which(rmd_content == "---")[2]
    rmd_yaml_content <- paste(rmd_content[(start_idx+1):(end_idx-1)], collapse = "\\n")
    rmd_params <- yaml::yaml.load(rmd_yaml_content)

    # Override the params
    rmd_params[['params']] <- modifyList(rmd_params[['params']], params)

    # Recursive function to add 'value' to list elements, except for top-level
    add_value_recursively <- function(lst, is_top_level = FALSE) {
        if (!is.list(lst)) {
            return(lst)
        }

        lst <- lapply(lst, add_value_recursively)
        if (!is_top_level) {
            lst <- list(value = lst)
        }
        return(lst)
    }

    # Reformat nested lists under 'params' to have a 'value' key recursively
    rmd_params[['params']] <- add_value_recursively(rmd_params[['params']], is_top_level = TRUE)

    # Convert back to YAML string
    updated_yaml_content <- as.character(yaml::as.yaml(rmd_params))

    # Remove the old YAML content
    rmd_content <- rmd_content[-((start_idx+1):(end_idx-1))]

    # Insert the updated YAML content at the right position
    rmd_content <- append(rmd_content, values = unlist(strsplit(updated_yaml_content, split = "\\n")), after = start_idx)

    writeLines(rmd_content, '${prefix}.parameterised.Rmd')

    # Render based on the updated file
    rmarkdown::render('${prefix}.parameterised.Rmd', output_file='${prefix}.html', envir = new.env())
    writeLines(capture.output(sessionInfo()), "session_info.log")
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmarkdown: \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch session_info.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmarkdown: \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))")
    END_VERSIONS
    """
}
