// IMPORTANT: Add this configuration to your modules.config

process {
    withName: ".*BAM_SPLIT_BY_REGION:SAMTOOLS_VIEW" {
        ext.args2 = {
            [
                "${meta.genomic_region}", // Specifies the chromosome name used for splitting the input bam file.
            ].join(' ').trim()
        }

        ext.prefix = { "${meta.id}_${meta.genomic_region}" }
    }

}
