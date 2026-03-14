/*
 * BIGSLICE_PREP_INPUT: Prepares antiSMASH outputs for BiG-SLiCE analysis
 * 
 * This process transforms antiSMASH output directories into the specific
 * directory structure and file format that BiG-SLiCE expects for clustering
 * biosynthetic gene clusters (BGCs).
 * 
 * BiG-SLiCE input structure:
 * input/
 * ├── datasets.tsv              # dataset configuration file
 * ├── <dataset_name>/           # GBK files organized by sample (each GBK contains a BGC region)
 * │   ├── sample1/
 * │   │   ├── region001.gbk     # BGC region 1 in GenBank format
 * │   │   └── region002.gbk     # BGC region 2 in GenBank format
 * │   └── sample2/
 * │       └── region001.gbk     # BGC region 1 in GenBank format
 * └── taxonomy/
 *     └── dataset_taxonomy.tsv  # taxonomic information (9-column format)
 */
process BIGSLICE_PREP_INPUT {
  label 'bigslice'
  tag "dataset=${params.bgc_bigslice_dataset_name}"

  conda "${moduleDir}/environment.yml"
container "${ workflow.containerEngine == 'singularity' && !task. ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0':
    'quay.io/biocontainers/bigslice:2.0.2--pyh8ed023e_0' }"

  input:
  // list of antiSMASH output directories (one per sample)
  path antismash_dirs

  output:
  // complete "input" folder structure for BiG-SLiCE (contains dataset/, taxonomy/, datasets.tsv)
  path "input", emit: input_dir
  path "versions.yml", emit: versions

  script:
  // prepare quoted directory list for bash for-loop processing
  def quoted = antismash_dirs.collect { "\"${it}\"" }.join(' ')
  def VERSION = '2.0.2'
  """
  set -euo pipefail

  ROOT="input"                              # BiG-SLiCE input root directory
  DS="${params.bgc_bigslice_dataset_name}"      # dataset name (e.g., 'antismash')
  OUT="\$ROOT/\$DS"                         # dataset-specific output directory
  TAXROOT="\$ROOT/taxonomy"                 # taxonomy directory

  rm -rf "\$ROOT"
  mkdir -p "\$OUT" "\$TAXROOT"

  for d in ${quoted}; do
    [ -d "\$d" ] || continue                # skip if directory doesn't exist
    sample=\$(basename "\$d")               # extract sample name from directory path
    mkdir -p "\$OUT/\$sample"               # create sample-specific subdirectory
    
    find -L "\$d" -type f \\( -name "*.region*.gbk" -o -name "*.gbk" \\) -print0 \
      | xargs -0 -I{} cp -f "{}" "\$OUT/\$sample/"
  done

  if [ -n "${params.bgc_bigslice_taxonomy ?: ''}" ]; then
    cp "${params.bgc_bigslice_taxonomy}" "\$TAXROOT/dataset_taxonomy.tsv"
  else
    printf "accession\\ttaxdomain\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\torganism\\n" > "\$TAXROOT/dataset_taxonomy.tsv"
    
    for d in "\$OUT"/*/; do
      [ -d "\$d" ] || continue
      acc=\$(basename "\$d")/                # sample accession (with trailing slash as per BiG-SLiCE format)
      printf "%s\\tUnknown\\tUnknown\\tUnknown\\tUnknown\\tUnknown\\tUnknown\\tUnknown\\tUnknown\\n" "\$acc" >> "\$TAXROOT/dataset_taxonomy.tsv"
    done
  fi

  {
    echo "# dataset_name\\tdataset_path\\ttaxonomy_path\\tdescription"
    printf "%s\\t%s\\t%s\\t%s\\n" "\$DS" "\$DS" "taxonomy/dataset_taxonomy.tsv" "antiSMASH \$DS"
  } > "\$ROOT/datasets.tsv"

    cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bigslice: \$(bigslice --version 2>&1 | grep -oP 'BiG-SLiCE \\K[0-9.]+' || echo "${VERSION}")
  END_VERSIONS
  """
}
