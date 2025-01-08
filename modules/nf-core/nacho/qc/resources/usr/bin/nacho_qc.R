#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(ggplot2)
library(fs)
library(NACHO)
library(readr)
library(tidyr)

# Commandline Argument parsing
option_list <- list(
    make_option(
        c("--input_rcc_path"),
        type = "character",
        default = "./" ,
        help = "Path to the folder that contains the RCC input file(s)",
        metavar = "character"),
    make_option(
        c("--input_samplesheet"),
        type = "character",
        default = NULL ,
        help = "Path to the sample sheet file",
        metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate mandatory arguments
if (is.null(opt$input_rcc_path)) {
    stop("Error: The --input_rcc_path parameter is mandatory and must be specified.")
}

if (is.null(opt$input_samplesheet)) {
    stop("Error: The --input_samplesheet parameter is mandatory and must be specified.")
}

input_rcc_path    <- opt$input_rcc_path
input_samplesheet <- opt$input_samplesheet

# Create filelist for NachoQC
list_of_rccs <- dir_ls(path = input_rcc_path, glob = "*.RCC")

# Core Code
nacho_data <- load_rcc(data_directory = input_rcc_path,
                        ssheet_csv = input_samplesheet,
                        id_colname = "RCC_FILE_NAME")

output_base <- "./"

# Write out HK genes detected and add to MultiQC report as custom content
line="#id: nf-core-nanostring-hk-genes
#section_name: 'Housekeeping Genes'
#description: 'The following Housekeeping Genes have been detected in the input RCC Files:'
#plot_type: 'html'
#section_href: 'https://github.com/nf-core/nanostring'
#data:
    "

write(line,file=paste0(output_base, "hk_detected_mqc.txt"),append=TRUE)
write(nacho_data$housekeeping_genes ,paste0(output_base,"hk_detected_mqc.txt"),append=TRUE)

# Add in all plots as MQC output for MultiQC
plot_bd <- autoplot(
    object = nacho_data,
    x = "BD",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="BD_mqc.png", plot_bd)

## Field of View (FoV) Imaging

plot_fov <- autoplot(
    object = nacho_data,
    x = "FoV",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="FOV_mqc.png", plot_fov)


## Positive Control Linearity

plot_posctrl_lin <- autoplot(
    object = nacho_data,
    x = "PCL",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)

ggsave(filename="Posctrl_linearity_mqc.png", plot_posctrl_lin)

## Limit of Detection

plot_lod <- autoplot(
    object = nacho_data,
    x = "LoD",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)

ggsave(filename="LOD_mqc.png", plot_lod)

## Positive Controls

plot_pos <- autoplot(
    object = nacho_data,
    x = "Positive",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="Pos_mqc.png", plot_pos)


## Negative Controls

plot_neg <- autoplot(
    object = nacho_data,
    x = "Negative",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="Neg_mqc.png", plot_neg)

## Housekeeping Genes

plot_hk <- autoplot(
    object = nacho_data,
    x = "Housekeeping",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="HK_mqc.png", plot_hk)

## Positive Controls vs Negative Controls

plot_pos_vs_neg <- autoplot(
    object = nacho_data,
    x = "PN",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="Pos_vs_neg_mqc.png", plot_pos_vs_neg)

## Average Counts vs. Binding Density

plot_avg_vs_bd <- autoplot(
    object = nacho_data,
    x = "ACBD",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="AVG_vs_BD_mqc.png", plot_avg_vs_bd)

## Average Counts vs. Median Counts

plot_avg_vs_med <- autoplot(
    object = nacho_data,
    x = "ACMC",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="AVG_vs_MED_mqc.png", plot_avg_vs_med)

## Principal Component 1 vs. 2

plot_pc12 <- autoplot(
    object = nacho_data,
    x = "PCA12",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="PCA1_vs_PCA2_mqc.png", plot_pc12)

## Principal Component i

plot_pcai <- autoplot(
    object = nacho_data,
    x = "PCAi",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="PCAi_mqc.png", plot_pcai)

## Principal Component planes
plot_pcap <- autoplot(
    object = nacho_data,
    x = "PCA",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="PCA_mqc.png", plot_pcap)

## Positive Factor vs. Negative Factor
plot_posf_vs_negf <- autoplot(
    object = nacho_data,
    x = "PFNF",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="POSF_vs_NEGF_mqc.png", plot_posf_vs_negf)

## Housekeeping Factor

plot_hkf <- autoplot(
    object = nacho_data,
    x = "HF",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="HKF_mqc.png", plot_hkf)

## Normalization Factors

plot_normf <- autoplot(
    object = nacho_data,
    x = "NORM",
    colour = "CartridgeID",
    size = 0.5,
    show_legend = TRUE
)
ggsave(filename="plot_normf_mqc.png", plot_normf)

# Create QC table for MultiQC Report
outliers_thresholds <- nacho_data[["outliers_thresholds"]]

qc_table <- nacho_data[["nacho"]] %>%
    select(c(RCC_FILE_NAME,BD,FoV,PCL,LoD,MC,MedC,Positive_factor,Negative_factor,House_factor)) %>%
    unique() %>%
    mutate("BD QC" = if_else(BD < outliers_thresholds[["BD"]][1] | BD > outliers_thresholds[["BD"]][2], "FAIL", "PASS"), .after = BD) %>%
    mutate("FoV QC" = if_else(FoV < outliers_thresholds[["FoV"]], "FAIL", "PASS"), .after = FoV) %>%
    mutate("PCL QC" = if_else(PCL < outliers_thresholds[["PCL"]], "FAIL", "PASS"), .after = PCL) %>%
    mutate("LoD QC" = if_else(LoD < outliers_thresholds[["LoD"]], "FAIL", "PASS"), .after = LoD) %>%
    mutate("PNF QC" = if_else(Positive_factor < outliers_thresholds[["Positive_factor"]][1] | Positive_factor > outliers_thresholds[["Positive_factor"]][2], "FAIL", "PASS"), .after = Positive_factor) %>%
    mutate("HKNF QC" = if_else(House_factor < outliers_thresholds[["House_factor"]][1] | House_factor > outliers_thresholds[["House_factor"]][2], "FAIL", "PASS"), .after = House_factor) %>%
    relocate(Negative_factor, .after = last_col()) %>%
    rename("Negative Factor" = Negative_factor) %>%
    rename("House Factor" = House_factor) %>%
    rename("Positive Factor" = Positive_factor) %>%
    rename("RCC_FILE" = RCC_FILE_NAME)

write_tsv(qc_table ,file=paste0(output_base,"normalized_qc_mqc.txt"))

# Render Standard Report for investigation in main MultiQC Report
render(nacho_data, output_dir = output_base, output_file = "NanoQC.html", show_outliers = FALSE)

# Render the same Report for standard investigation, but not for MultiQC Report
render(nacho_data, output_dir = output_base, output_file = "NanoQC_with_outliers.html", show_outliers = TRUE)
