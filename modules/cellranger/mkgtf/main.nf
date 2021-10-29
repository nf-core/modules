// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CELLRANGER_MKGTF {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "nfcore/cellranger:6.0.2"

    input:
    path gtf

    output:
    path "versions.yml"   , emit: versions
    path "*.filtered.gtf" , emit: gtf

    script:
    """
    cellranger mkgtf $gtf ${gtf.baseName}.filtered.gtf \
                    --attribute=gene_biotype:protein_coding \
                    --attribute=gene_biotype:lincRNA \
                    --attribute=gene_biotype:antisense \
                    --attribute=gene_biotype:IG_LV_gene \
                    --attribute=gene_biotype:IG_V_gene \
                    --attribute=gene_biotype:IG_V_pseudogene \
                    --attribute=gene_biotype:IG_D_gene \
                    --attribute=gene_biotype:IG_J_gene \
                    --attribute=gene_biotype:IG_J_pseudogene \
                    --attribute=gene_biotype:IG_C_gene \
                    --attribute=gene_biotype:IG_C_pseudogene \
                    --attribute=gene_biotype:TR_V_gene \
                    --attribute=gene_biotype:TR_V_pseudogene \
                    --attribute=gene_biotype:TR_D_gene \
                    --attribute=gene_biotype:TR_J_gene \
                    --attribute=gene_biotype:TR_J_pseudogene \
                    --attribute=gene_biotype:TR_C_gene

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
