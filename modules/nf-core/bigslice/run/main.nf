/*
 * BIGSLICE_RUN: Executes BiG-SLiCE clustering analysis on prepared BGC data
 * 
 * BiG-SLiCE (Biosynthetic Gene cluster - Super Linear Clustering Engine) is a tool
 * for rapid clustering and analysis of biosynthetic gene clusters from genomic data.
 * 
 * This process:
 * 1. Takes the structured input directory created by BIGSLICE_PREP_INPUT
 * 2. Uses pre-trained machine learning models to analyze BGC features
 * 3. Performs hierarchical clustering of BGCs based on genetic similarity
 * 4. Generates interactive web application and detailed clustering results
 * 
 * Output structure:
 * output/
 * ├── app/                           # interactive web application for visualization
 * ├── result/                        # clustering results and analysis data
 * ├── LICENSE.txt                    # BiG-SLiCE license information
 * ├── requirements.txt               # Python dependencies for the web app
 * └── start_server.sh               # script to launch the interactive web server
 */
process BIGSLICE_RUN {
  label 'bigslice'

  conda "${moduleDir}/enviroment.yml"
container "${ workflow.containerEngine == 'singularity' && !task. ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0':
    'quay.io/biocontainers/bigslice:2.0.2--pyh8ed023e_0' }"

  input:
  path input_dir    // structured input directory from BIGSLICE_PREP_INPUT (contains datasets.tsv, BGC files, taxonomy)
  path models_dir   // pre-trained BiG-SLiCE models directory (e.g., bigslice-models.2022-11-30)

  output:
  path "output", emit: outdir   // complete BiG-SLiCE output directory with clustering results and reports
  path "versions.yml", emit: versions

  script:
  def VERSION = '2.0.2'
  """
  set -euo pipefail
  
  # clean any existing output directory to ensure fresh results
  rm -rf output 2>/dev/null || true

  # execute BiG-SLiCE clustering analysis
  # -i: input directory containing prepared BGC data and configuration
  # --program_db_folder: directory with pre-trained machine learning models
  # output: destination directory for all results and reports
  bigslice \
    -i "${input_dir}" \
    --program_db_folder "${models_dir}" \
    output
      
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bigslice: \$(bigslice --version 2>&1 | grep -oP 'BiG-SLiCE \\K[0-9.]+' || echo "${VERSION}")
  END_VERSIONS
  """
}
