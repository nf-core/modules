process LAST_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0d/0d27a2649f1291ff817dc8f73357ffac206424cd972d3855421e4258acc600f7/data'
        : 'community.wave.seqera.io/library/last:1611--e1193b3871fa0975'}"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    tuple val(meta), path("*.tsv")   , emit: multiqc
    // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    set -o pipefail

    # LAST reports genome sizes and sequence number at the beginning and end of the MAF files it outputs.
    function get_genome_stats() { awk '
        BEGIN { OFS = "\\t" }
        function extract(line, regex,   s) {
            if (match(line, regex)) {
                s = substr(line, RSTART, RLENGTH)
                sub(/^[^=]*=/, "", s)   # strip up to "="
                return s
            }
            return ""
        }
        /^# Reference sequences=/ {
            ref_seq     = extract(\$0, "^# Reference sequences=[0-9]+")
            ref_letters = extract(\$0,                "letters=[0-9]+")
        }
        /^# Query sequences=/ {
            qry_seq     = extract(\$0, "^# Query sequences=[0-9]+")
            qry_letters = extract(\$0,            "letters=[0-9]+")
        }
        END {
            print "TargetSequences", "TargetLength", "QuerySequences", "QueryLength"    # Header for MultiQC
            print ref_seq+0,         ref_letters+0,  qry_seq+0,         qry_letters+0   # Data in TSV format
        }'
    }

    # The MAF files do not report number of matches directly, but we can compute it easily via conversion to PSL format and then extract them.
    function calculate_psl_metrics() {
        awk 'BEGIN {
            OFS="\\t";
            FS="\\t";  # Set field separator as tab
            totalMatches = 0;
            totalAlignmentLength = 0;
            totalAlignedBases = 0;
        }
        {
            totalMatches         += \$1 +       \$3            ;  # Sum matches          and repMatches
            totalAlignmentLength += \$1 + \$2 + \$3 + \$6 + \$8;  # Sum matches, misMatches, repMatches, qBaseInsert, and tBaseInsert
            totalAlignedBases    += \$1 + \$2 + \$3            ;  # Sum matches, misMatches, repMatches
        }
        END {
            percentIdentity       = (totalAlignmentLength > 0) ? (totalMatches / totalAlignmentLength * 100) : 0;
            percentIdentityNoGaps = (totalAlignmentLength > 0) ? (totalMatches / totalAlignedBases    * 100) : 0;
            print "Sample",  "TotalAlignmentLength", "PercentIdentity", "PercentIdentityNoGaps";  # Header for MultiQC
            print "$meta.id", totalAlignmentLength,   percentIdentity,   percentIdentityNoGaps;   # Data in TSV format
        }'
    }

    # The MAF files can be really big, so we stream them in the awk functions and gzip instead of reading them each time.
    zcat < $maf |
        last-split $args |
        tee >(get_genome_stats > ${prefix}.genomestats.txt) |
        tee >(gzip --no-name   > ${prefix}.maf.gz) |
        maf-convert psl |
        calculate_psl_metrics  > ${prefix}.alignmentstats.txt

    # Combine the two stats file into one for MultiQC.
    paste ${prefix}.alignmentstats.txt ${prefix}.genomestats.txt > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo stub | gzip --no-name > ${prefix}.maf.gz
    touch ${prefix}.tsv
    """
}
