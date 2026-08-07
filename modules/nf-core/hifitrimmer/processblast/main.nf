process HIFITRIMMER_PROCESSBLAST {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b30f66d13fc3da362a0a9abb9365a6b5e2b90e6aff316c59ae58471c707bdd6/data' :
      'community.wave.seqera.io/library/hifi_trimmer:4.0.0--41f93366a6098f77' }"

   input:
   tuple val(meta), path(blast)
   tuple val(meta2), path(yaml)


   output:

   tuple val(meta), path("*.bed.gz")      , emit: bed
   tuple val(meta), path("*.summary.json"), emit: summary
   tuple val(meta), path("*.hits")        , emit: hits, optional: true
   tuple val("${task.process}"), val('hifi-trimmer'), eval("hifi-trimmer --version | cut -d' ' -f3"), emit: versions_hifitrimmer, topic: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ? task.ext.args : ''
   """
   hifi-trimmer process-blast \\
      -t ${task.cpus} \\
      --prefix ${prefix} \\
      ${args} \\
      ${blast} \\
      ${yaml}
   """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   # Create deterministic gzip output (fixed mtime) so stub md5 is stable across runs.
   python3 -c "import gzip; open('${prefix}.bed.gz','wb').write(gzip.compress(b'stub\\n', mtime=0))"
   touch ${prefix}.summary.json
   touch ${prefix}.hits
   """
}
