// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
include { dump_params_yml } from "./parametrize"

params.options = [:]
options        = initOptions(params.options)
params.parametrize = true
params.implicit_params = true
params.meta_params = true

process RMARKDOWNNOTEBOOK {
    // tag { meta.id }
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    //NB: You likely want to override this with a container containing all required
    //dependencies for you analysis. The container at least needs to contain the
    //yaml and rmarkdown R packages.
    conda (params.enable_conda ? "r-base=4.1.0 r-rmarkdown=2.9 r-yaml=2.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        // TODO replace with biocontainer
        container "https://depot.galaxyproject.org/singularity/TODO"
    } else {
        container "quay.io/biocontainers/TODO"
    }

    input:
    tuple val(meta), path(notebook)
    val(parameters)
    path(input_files)

    output:
    tuple val(meta), path("*.html"), emit: report
    path("artifacts/*"), emit: artifacts, optional: true
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def params_cmd = ""
    def render_cmd = ""
    if (params.parametrize) {
        nb_params = [:]
        if (params.implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "./"
        }
        if (params.meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += parameters
        params_cmd = dump_params_yml(nb_params)
        render_cmd = """
            params = yaml::read_yaml('.params.yml')
            rmarkdown::render('${prefix}.Rmd', params=params, envir=new.env())
        """
    } else {
        render_cmd = "rmarkdown::render('${prefix}.Rmd')"
    }

    params_cmd +
    """
    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"

    # Work around  https://github.com/rstudio/rmarkdown/issues/1508
    # If the symbolic link is not replaced by a physical file
    # output- and temporary files will be written to the original directory.
    mv "${notebook}" "${notebook}.orig"
    cp -L "${notebook}.orig" "${prefix}.Rmd"

    # Render notebook
    Rscript - <<EOF
    ${render_cmd}
    EOF

    echo \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))") > ${software}.version.txt
    """.stripIndent()
}
