#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Local default params
params.internal_outdir = 'results'
params.internal_process_name = 'dedup'

// Check for internal parameter overrides
nfUtils.check_internal_overrides(module_name, params)

/*-------------------------------------------------> DEDUP OPTIONS <-----------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------
EXTRACT BARCODE OPTION -> --extract-umi-method
-------------------------------------------------------------------------------------------------------------------------------*/

//Default extract barcode -> Barcodes are contained at the end of the read separated as specified with --umi-separator option
//Enter the character that corresponds to the UMI separator
params.internal_umi_separator = ':'

/*-----------------------------------------------------------------------------------------------------------------------------*/

//Barcodes contained in a tag(s) -> --extract-umi-method=tag
params.internal_extract_method_tag = false

// IF params.internal_extract_method_tag = true -> choose between the 6 options below
//Umi tag options
//--umi-tag=[TAG], --umi-tag-split=[SPLIT], --umi-tag-delimiter=[DELIMITER]
//Select tag, split ot delimiter
params.internal_umi_tag_active = false
params.internal_umi_tag_split_active = false
params.internal_umi_tag_delim_active = false
//If one set to true, insert tag, split or delimiter in the param below
params.internal_umi_tag = ''


//Cell tag options
//--umi-tag=[TAG], --umi-tag-split=[SPLIT], --umi-tag-delimiter=[DELIMITER]
//Select tag, split ot delimiter
params.internal_cell_tag_active = false
params.internal_cell_tag_split_active = false
params.internal_cell_tag_delim_active = false
//If one set to true, insert tag, split or delimiter in the param below
params.internal_cell_tag = ''

/*-----------------------------------------------------------------------------------------------------------------------------
OUTPUT OPTION -> --stdout or -S
-------------------------------------------------------------------------------------------------------------------------------*/
//Insert output file name
params.internal_output_file_name = ''
//Activate this option to use sample_id in the output file names -> e.g. sample1.dedup.bam
params.internal_output_sampleid = true

/*-----------------------------------------------------------------------------------------------------------------------------
OUTPUT STATS OPTION -> --output-stats=[PREFIX]
-------------------------------------------------------------------------------------------------------------------------------*/
//Enter the prefix 
params.internal_output_stats = ''
//Activate this option to use sample_id as the prefix -> e.g. sample1_edit_distance
params.internal_output_stats_sampleid = true

/*-----------------------------------------------------------------------------------------------------------------------------
GROUPING METHOD OPTION -> --method=[method]
-------------------------------------------------------------------------------------------------------------------------------*/
//What method to use to identify group of reads with the same (or similar) UMI(s)
//Default method is directional
//Choose a grouping method: unique, percentile, cluster or adjacency

//Reads group share the exact same UMI
params.internal_grouping_unique = false

//Reads group share the exact same UMI. UMIs with counts < 1% of the median counts for UMIs at the same position are ignored.
params.internal_grouping_percentile = false

//Identify clusters of connected UMIs (based on hamming distance threshold). Each network is a read group
params.internal_grouping_cluster = false

//Cluster UMIs as above. For each cluster, select the node (UMI) with the highest counts.
params.internal_grouping_adjacency = false

/*-----------------------------------------------------------------------------------------------------------------------------*/
//Additional grouping method options

//--edit-threshold-distance
//For the adjacency and cluster methods, the threshold for the edit distance to connect two UMIs in the network can be increased. The default value of 1 works best unless the UMI is very long (>14bp).
params.internal_edit_threshold_distance = ''

//--spliced-is-unique
//Causes two reads that start in the same position on the same strand and having the same UMI to be considered unique if one is spliced and the other is not.
params.internal_unique_spliced = false

//--soft-clip-threshold
//By setting this option, you can treat reads with at least this many bases soft-clipped at the 3’ end as spliced. Default=4
params.internal_soft_clip_threshold = ''

//--read-length
//Use the read length as a criteria when deduping, for e.g sRNA-Seq.
params.internal_read_length = true

/*-----------------------------------------------------------------------------------------------------------------------------
SINGLE-CELL RNA-SEQ OPTIONS 
-------------------------------------------------------------------------------------------------------------------------------*/
//--per-gene
//Reads will be grouped together if they have the same gene.
//MUST SUPPLY --per-contig OR --gene-tag
params.internal_per_gene = false

//--gene-tag
//Deduplicate per gene. The gene information is encoded in the bam read tag specified
// USE WITH --per-gene
params.internal_gene_tag = false

//--assigned-status-tag
//BAM tag which describes whether a read is assigned to a gene. Defaults to the same value as given for --gene-tag
params.internal_assigned_status_tag = ''

//--skip-tags-regex
//Use in conjunction with the --assigned-status-tag option to skip any reads where the tag matches this regex. 
//Default ("^[__|Unassigned]") matches anything which starts with “__” or “Unassigned”:
params.internal_skip_tags_regex = ''

//--per-contig
//Deduplicate per contig (field 3 in BAM; RNAME). All reads with the same contig will be considered to have the same alignment position. 
params.internal_per_contig = false

//--gene-transcript-map
//File mapping genes to transcripts (tab separated)
params.internal_gene_transcript_map = false

//--per-cell
//Reads will only be grouped together if they have the same cell barcode. Can be combined with --per-gene.
params.internal_per_cell = false

/*-----------------------------------------------------------------------------------------------------------------------------
SAM/BAM OPTIONS -> --method=[method]
-------------------------------------------------------------------------------------------------------------------------------*/

//--mapping-quality
//Minimium mapping quality (MAPQ) for a read to be retained. Default is 0.
params.internal_mapping_quality = 0

//--unmapped-reads -> How should unmapped reads be handled
// Activate one of those three options: discard (default), use, output
params.internal_unmapped_reads_discard = true 
params.internal_unmapped_reads_use = false 
params.internal_unmapped_reads_output = false 


//--chimeric-pairs -> How should chimeric pairs be handled
// Activate one of those three options: discard, use (default), output
params.internal_chimeric_pairs_discard = false
params.internal_chimeric_pairs_use = false
params.internal_chimeric_pairs_output = false


//--unpaired-reads -> How should unpaired reads be handled
// Activate one of those three options: discard, use (default), output
params.internal_unpaired_reads_discard = false 
params.internal_unpaired_reads_use = false 
params.internal_unpaired_reads_output = false 

//--ignore-umi
//Ignore the UMI and group reads using mapping coordinates only
//CANNOT BE USED WITH --output-stats
params.internal_ignore_umi = false

//--subset
//Only consider a fraction of the reads, chosen at random. This is useful for doing saturation analyses.
params.internal_subset = ''

//--chrom
//Only consider a single chromosome. This is useful for debugging/testing purposes
params.internal_chrom = false

// Input/output BAM options

// Parameters for input and output in SAM format instead of BAM
params.internal_in_sam = false 
params.internal_out_sam = false

//BAM is paired end - output both read pairs. This will also force the use of the template length to determine reads with the same mapping coordinates.
//--paired
params.internal_paired_end = false

/*-----------------------------------------------------------------------------------------------------------------------------
GROUP/DEDUP OPTIONS 
-------------------------------------------------------------------------------------------------------------------------------*/

//--no-sort-output
params.internal_no_sort_output = false

//--buffer-whole-contig
params.internal_buffer_whole_contig = false

/*-----------------------------------------------------------------------------------------------------------------------------*/




// dedup reusable component
process dedup {
    publishDir "umi_tools/dedup/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-umitools:latest'

    input:
      tuple val(sample_id), path(bai), path(bam)
       
    output:
      tuple val(sample_id), path(bam), emit: dedupBam

    script:

    //Initializing dedup arguments
    dedup_pre_args = "umi_tools dedup "
    dedup_post_args = ''

    //Method used to extract barcodes
    if (params.internal_umi_separator != ''){
      dedup_pre_args += "--umi-separator=\"$params.internal_umi_separator\" "
    }

    //Input args 
    dedup_pre_args += "-I "

    //Output_stats option
    if (params.internal_output_stats_sampleid){
      dedup_post_args += "--output-stats=$sample_id "
    } else {
      if (params.internal_output_stats != ''){
      dedup_post_args += "--output-stats=$params.internal_output_stats "
      }
    }

    //Output/ stdout option 
    if (params.internal_output_sampleid){
      dedup_post_args += "-S ${sample_id}.dedup.bam "
    }else {
      if (params.internal_output_file_name != ''){
        dedup_post_args += "-S $params.internal_output_file_name "
      }
    }

    //Grouping method option 
    if (params.internal_grouping_unique){
      dedup_post_args += "--method=unique "
    }
    if (params.internal_grouping_percentile){
      dedup_post_args += "--method=percentile "
    }
    if (params.internal_grouping_cluster){
      dedup_post_args += "--method=cluster "
    }
    if (params.internal_grouping_adjacency){
      dedup_post_args += "--method=adjacency "
    }

    //Additional grouping method options 
    if (params.internal_edit_threshold_distance != ''){
      dedup_post_args += "--edit-threshold-distance=$params.internal_edit_threshold_distance "
    }
    if (params.internal_unique_spliced){
      dedup_post_args += "--spliced-is-unique "
    }
    if (params.internal_soft_clip_threshold != ''){
      dedup_post_args += "--soft-clip-threshold=$params.internal_soft_clip_threshold "
    }
    if (params.internal_read_length){
      dedup_post_args += "--read-length "
    }

    //Single-cell RNA-seq options
    if (params.internal_per_gene){
      dedup_post_args += "--per-gene "
    }
    if (params.internal_gene_tag){
      dedup_post_args += "--gene-tag "
    }
    if (params.internal_assigned_status_tag != ''){
      dedup_post_args += "--assigned-status-tag=$params.internal_assigned_status_tag "
    }
    if (params.internal_skip_tags_regex != ''){
      dedup_post_args += "--skip-tags-regex=$params.internal_skip_tags_regex "
    }
    if (params.internal_per_contig){
      dedup_post_args += "--per-contig "
    }
    if (params.internal_gene_transcript_map){
      dedup_post_args += "--gene-transcript-map "
    }
    if (params.internal_per_cell){
      dedup_post_args += "--per-cell "
    }

    // SAM/BAM options
    // Unmapped reads options
    if (params.internal_unmapped_reads_discard){
      dedup_post_args += "--unmapped-reads "
    } 
    if (params.internal_unmapped_reads_use){
      dedup_post_args += "--unmapped-reads=use "
    }
    if (params.internal_unmapped_reads_output){
      dedup_post_args += "--unmapped-reads=output "
    }

    // Chimeric pairs options
    if (params.internal_chimeric_pairs_discard){
      dedup_post_args += "--chimeric-pairs=discard "
    }
    if (params.internal_chimeric_pairs_use){
      dedup_post_args += "--chimeric-pairs "
    }
    if (params.internal_chimeric_pairs_output){
      dedup_post_args += "--chimeric-pairs=output "
    }

    // Unpaired reads options
    if (params.internal_unpaired_reads_discard){
      dedup_post_args += "--unpaired-reads=discard "
    } 
    if (params.internal_unpaired_reads_use){
      dedup_post_args += "--unpaired-reads "
    }
    if (params.internal_unpaired_reads_output){
      dedup_post_args += "--unpaired-reads=output "
    }

    // Additional SAM/BAM options
    if (params.internal_mapping_quality != 0){
      dedup_post_args += "--mapping-quality=$params.internal_mapping_quality "
    }
    if (params.internal_ignore_umi){
      dedup_post_args += "--ignore-umi "
    }
    if (params.internal_subset != ''){
      dedup_post_args += "--subset=$params.internal_subset "
    }
    if (params.internal_chrom){
      dedup_post_args += "--chrom "
    }

    // Input/output BAM options
    if (params.internal_in_sam){
      dedup_post_args += "--in-sam "
    }
    if (params.internal_out_sam){
      dedup_post_args += "--out-sam "
    }
    if (params.internal_paired_end){
      dedup_post_args += "--paired "
    }

    //Group dedup options
    if (params.internal_no_sort_output){
      dedup_post_args += "--no-sort-output "
    }
    if (params.internal_buffer_whole_contig){
      dedup_post_args += "--buffer-whole-contig "
    }

    // Displays the umi_tools command line to check for mistakes
    println dedup_pre_args
    println dedup_post_args

    """
    $dedup_pre_args $bam $dedup_post_args 
    """
}