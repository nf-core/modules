process CELLRANGER_MULTI {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/nfcore/cellranger:7.1.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_MULTI module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
	val meta
    tuple val(meta_gex)    , path (gex_fastqs   , stageAs: "fastqs/gex/*")
    tuple val(meta_vdj)    , path (vdj_fastqs   , stageAs: "fastqs/vdj/*")
    tuple val(meta_ab)     , path (ab_fastqs    , stageAs: "fastqs/ab/*")
    tuple val(meta_beam)   , path (beam_fastqs  , stageAs: "fastqs/beam/*")
    tuple val(meta_cmo)    , path (cmo_fastqs   , stageAs: "fastqs/cmo/*")
    tuple val(meta_crispr) , path (crispr_fastqs, stageAs: "fastqs/crispr/*")
    path gex_reference     , stageAs: "references/gex/*"
    path gex_frna_probeset , stageAs: "references/gex/probeset/*"
    path gex_targetpanel   , stageAs: "references/gex/targetpanel/*"
    path vdj_reference     , stageAs: "references/vdj/*"
    path vdj_primer_index  , stageAs: "references/vdj/primers/*"
    path fb_reference      , stageAs: "references/fb/*"
    path beam_panel        , stageAs: "references/beam/panel/*"
    path cmo_reference     , stageAs: "references/cmo/*"
    path cmo_barcodes      , stageAs: "references/cmo/barcodes/*"
    path frna_sampleinfo   , stageAs: "references/frna/*"

    output:
    tuple val(meta), path("cellranger_multi_config.csv"), emit: config
    tuple val(meta), path("${meta.id}/outs/*")          , emit: outs
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

	// if references + FASTQ are empty, then don't run corresponding analyses
	// get names of references, if they exist
    // empty reference channels stage as "references" 
    // empty FASTQ channels stage as "fastqs"
    // empty files stage as the file name, we check against 'EMPTY'
	def gex_reference_name     = gex_reference.getName() != 'references'     ? gex_reference.getName() : ''
	def gex_frna_probeset_name = gex_frna_probeset.getBaseName() != 'EMPTY'  ? gex_frna_probeset.getName() : ''
	def gex_targetpanel_name   = gex_targetpanel.getBaseName() != 'EMPTY'    ? gex_targetpanel.getName() : ''
	def vdj_reference_name     = vdj_reference.getName() != 'references'     ? vdj_reference.getName() : ''

	def include_gex  = gex_fastqs.first().getName() != 'fastqs' && gex_reference ? "[gene-expression]\nreference,references/gex/${gex_reference_name}" : ''
	def include_vdj  = vdj_fastqs.first().getName() != 'fastqs' && vdj_reference ? "[vdj]\nreference,references/vdj/${vdj_reference_name}" : ''
	def include_beam = beam_fastqs.first().getName() != 'fastqs' && beam_panel   ? '[antigen-specificity]' : ''
	def include_cmo  = cmo_fastqs.first().getName() != 'fastqs' && cmo_barcodes  ? '[samples]' : ''
	def include_fb   = fb_reference.first().getName() != 'references'            ? '[feature]' : ''
	def include_frna = gex_frna_probeset_name && frna_sampleinfo                 ? '[samples]' : ''

    println("fb_reference: ${fb_reference.first().getName()}")
    println("beam_fastqs: ${beam_fastqs.first().getName()}")
    println("cmo_fastqs: ${cmo_fastqs.first().getName()}")
    println("frna_fastqs: ${gex_frna_probeset.first().getName()}")
    println("include_beam: $include_beam")
    println("include_cmo: $include_cmo")
    println("include_fb: $include_fb")
    println("include_frna: $include_frna")
    println("gex_frna_probeset_name: actual ${gex_frna_probeset.getBaseName()}, expected $gex_frna_probeset_name")
    println("gex_targetpanel_name: actual ${gex_targetpanel.getBaseName()}, expected $gex_targetpanel_name")

	// targeted GEX panel goes under GEX section, not its own
	//def target_panel   = gex_targetpanel_name != 'references' ? "target-panel,\$PWD/references/$gex_targetpanel_name" : '' 
	def target_panel   = gex_targetpanel_name != '' ? "target-panel,\$PWD/references/gex/targetpanel/$gex_targetpanel_name" : '' 
    println("$target_panel")
    println("${gex_targetpanel_name != '' ? 'check target panel name!!!' : 'ok'}")

	// fixed RNA reference (not sample info!) also goes under GEX section
	//def frna_probeset  = include_frna && gex_frna_probeset_name ? "probe-set,\$PWD/references/gex_probeset" : ''
	def frna_probeset  = include_frna && gex_frna_probeset_name ? "probe-set,\$PWD/references/gex/probeset/$gex_frna_probeset_name" : ''
    println("$frna_probeset")
    println("${gex_frna_probeset_name != '' ? 'check frna probeset name!!!' : 'ok'}")

	// VDJ inner primer set
	def primer_index   = vdj_primer_index.getBaseName() != 'EMPTY' ? "inner-enrichment-primers,\$PWD/references/primers/${vdj_primer_index.getName()}" : '' 
    println("$primer_index")
    println("vdj_primer_index: actual $primer_index, expected ${vdj_primer_index.getBaseName()}")
    println("${primer_index != '' ? 'check primer index name!!!' : 'ok'}")

	// collect options for each section
	// these are pulled from the meta maps
	//def gex_options    = include_gex && meta_gex?.options ? meta_gex.options.inject('') { key, value -> key ? "$key,$value\n" : '' } : ''
	//def gex_options    = include_gex && meta_gex?.options ? meta_gex.options.inject('') { key, value -> value ? "${value.replaceFirst(/'='/, ",")}\n" : '' } : ''
	def gex_options    = include_gex && meta_gex?.options ? meta_gex.options.collect{ key, value -> /$key,$value/ }.join("\n") : ''
	def vdj_options    = include_vdj && meta_vdj?.options ? meta_vdj.options.inject('') { key, value -> key ? "$key,$value\n" : '' } : ''
	def ab_options     = include_fb && meta_ab?.options ? meta_ab.options.inject('') { key, value -> key ? "$key,$value\n" : '' } : ''
	def beam_options   = include_beam && meta_beam?.options ? meta_beam.options.inject('') { key, value -> key ? "$key,$value\n" : '' } : ''
	def cmo_options    = include_cmo && meta_cmo?.options ? meta_cmo.options.inject('') { key, value -> key ? "$key,$value\n" : '' } : ''
	def crispr_options = include_fb && meta_crispr?.options ? meta_crispr.options.inject('') { key, value -> key ? "$key,$value\n" : '' } : ''

	// the feature barcodes section get options for either CRISPR or antibody capture assays
	def fb_options     = (ab_options?.trim()) ? ab_options : crispr_options 

	// pull CSV text from these reference panels
	// these references get appended directly to config file
    println("gex_options: $gex_options")
    println("meta_gex: $meta_gex")
    println("meta_gex?.options: ${meta_gex?.options}")
    println("${meta_gex.options["expect-cells"]}")
    println("${meta_gex.options.chemistry}")
    println("vdj_options: $vdj_options")
    println("ab_options: $ab_options")
    println("beam_options: $beam_options")
    println("cmo_options: $cmo_options")
    println("crispr_options: $crispr_options")
    println("$beam_panel")
    println("$cmo_barcodes")
    println("$frna_sampleinfo")
	def beam_csv_text = include_beam && beam_panel.size() > 0      ? beam_panel.text : ''
	def cmo_csv_text  = include_cmo  && cmo_barcodes.size() > 0    ? cmo_barcodes.text : ''
	def frna_csv_text = include_frna && frna_sampleinfo.size() > 0 ? frna_sampleinfo.text : ''
    println("$beam_csv_text")
    println("$cmo_csv_text")
    println("$frna_csv_text")

	// point config to FASTQs
	def fastq_gex      = include_gex                  ? "${meta_gex.id},\$PWD/fastqs/gex/,,Gene Expression" : ''
	def fastq_vdj      = include_vdj                  ? "${meta_vdj.id},\$PWD/fastqs/vdj,,VDJ" : ''
	def fastq_antibody = include_fb && ab_options     ? "${meta_ab.id},\$PWD/fastqs/ab,,Antibody Capture" : ''
	def fastq_beam     = include_beam                 ? "${meta_beam.id},\$PWD/fastqs/beam,,Antigen Capture" : ''
	def fastq_crispr   = include_fb && crispr_options ? "${meta_crispr.id},\$PWD/fastqs/crispr,,CRISPR Guide Capture" : ''
	def fastq_cmo      = include_cmo                  ? "${meta_cmo.id},\$PWD/fastqs/cmo,,Multiplexing Capture" : ''

    println("fastq_gex: $fastq_gex")
    println("fastq_vdj: $fastq_vdj")
    println("fastq_antibody: $fastq_antibody")
    println("fastq_beam: $fastq_beam")
    println("fastq_crispr: $fastq_crispr")
    println("fastq_cmo: $fastq_cmo")

    """
    # stage cellranger multi config
	config="cellranger_multi_config.csv"
	touch \$config

	# add GEX if requested
	if [ "$include_gex" ]; then
		echo "\n$include_gex\n$gex_options" >> \$config

		# add CMO if included
		if [ "$cmo_options" ]; then echo "$cmo_options" >> \$config; fi
		
		# add fRNA if included
		if [ "$frna_probeset" ]; then echo "$frna_probeset" >> \$config; fi

		# add tGEX if included
		if [ "$target_panel" ]; then echo "$target_panel" >> \$config; fi
	fi

	# add VDJ if requested
	if [ "$include_vdj" ]; then
		echo "\n$include_vdj\n$vdj_options" >> \$config

		# add inner primer set if included
		if [ "$primer_index" ]; then echo "$primer_index" >> \$config; fi
	fi

	# add feature barcodes if requested
	if [ "$include_fb" ]; then
		echo "\n$include_fb\n$fb_options" >> \$config
	fi

	# add antigen capture if requested
	if [ "$include_beam" ]; then
		echo "\n$include_beam\n$beam_csv_text" >> \$config
	fi

	# add cell multiplexing if requested
	if [ "$include_cmo" ]; then
		echo "\n$include_cmo\n$cmo_csv_text" >> \$config
	fi

	# add fixed RNA profiling if requested
	if [ "$include_frna" ]; then
		echo "\n$include_frna\n$frna_csv_text" >> \$config
	fi

	# list FASTQ directories
	echo "[libraries]" >> \$config
    #echo "fastq_id,fastqs,lanes,feature_types,subsample_rate" >> \$config
    echo "fastq_id,fastqs,lanes,feature_types" >> \$config
	echo "$fastq_gex" >> \$config
	echo "$fastq_vdj" >> \$config
	echo "$fastq_antibody" >> \$config
	echo "$fastq_beam" >> \$config
	echo "$fastq_crispr" >> \$config
	echo "$fastq_cmo" >> \$config

    # output existing config to stdout
    echo "--- beginning of cellranger multi config ---"
    cat \$config
    echo "--- end of cellranger multi config ---"

	# run cellranger using config written above
    cellranger \\
        multi \\
        --id="${meta.id}" \\
        --csv=\$config \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}/outs/"
    touch ${meta.id}/outs/fake_file.txt
	touch cellranger_multi_config.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
