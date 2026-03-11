#!/usr/bin/env Rscript

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    args_vals = lapply(args_list, function(x) {
        x_splt = strsplit(x, split=":")[[1]]
        c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
    })

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]


# Load packages

suppressPackageStartupMessages({
library(ctree)
library(mobster)
library(VIBER)
library(cli)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
})

outdir = ""

# Auxiliary functions

add_dummy_driver = function(input_table, variant_colname, is_driver_colname) {
    add_driver = FALSE
    idx = dplyr::case_when(
        "cluster" %in% colnames(input_table) ~ which(input_table[["cluster"]] != "Tail")[1],
        .default = 1
    )
    if (!variant_colname %in% colnames(input_table) | !is_driver_colname %in% colnames(input_table)) {
        input_table = input_table %>% dplyr::mutate(!!is_driver_colname:=FALSE, !!variant_colname:=NA)
        add_driver = TRUE
    } else if (all(input_table[[is_driver_colname]]==FALSE)) {
        add_driver = TRUE
    }
    if (add_driver) {
        input_table[idx, is_driver_colname] = TRUE
        input_table[idx, variant_colname] = ""
    }
    return(input_table)
}

initialize_ctree_obj_pyclone = function(ctree_input) {
    driver_cluster = unique(ctree_input[which(ctree_input["is.driver"]==TRUE),c("cluster")])
    # the CCF table must report CCF values for each cluster and sample
    # cluster | nMuts | is.driver | is.clonal | sample1 | sample2 | ...
    CCF_table = ctree_input %>%
        dplyr::select(sample_id, cluster, nMuts, is.driver, is.clonal, CCF) %>%
        dplyr::mutate(is.driver=ifelse(is.driver=="", FALSE, TRUE)) %>%
        dplyr::filter(cluster!="Tail") %>%
        dplyr::group_by(cluster) %>%
        dplyr::mutate(is.driver=any(is.driver)) %>%
        dplyr::filter(any(CCF>0)) %>%
        dplyr::ungroup() %>%
        unique() %>%
        tidyr::pivot_wider(names_from="sample_id", values_from="CCF", values_fill=0)

    # the driver table must contain patient and variant IDs and report clonality and driver status
    # patientID | variantID | is.driver | is.clonal | cluster | sample1 | sample2 | ...
    drivers_table = ctree_input %>%
        dplyr::filter(cluster %in% CCF_table[["cluster"]]) %>%
        dplyr::mutate(is.driver=as.logical(is.driver)) %>%
        dplyr::select(patientID, sample_id, variantID, cluster, is.driver, is.clonal, CCF) %>%
        dplyr::filter(is.driver==TRUE) %>%
        dplyr::mutate(variantID=replace(variantID, is.na(variantID), "")) %>%
        tidyr::pivot_wider(names_from="sample_id", values_from="CCF", values_fill=0)

    samples = unique(ctree_input[["sample_id"]])  # if multisample, this is a list
    patient = unique(ctree_input[["patientID"]])

    CCF_table = add_dummy_driver(CCF_table, variant_colname="variantID", is_driver_colname="is.driver") %>%
        dplyr::mutate(cluster=as.character(cluster))

    if (nrow(drivers_table)==0) {
        drivers_table = CCF_table %>%
            dplyr::filter(is.driver) %>%
            dplyr::select(-nMuts) %>%
            dplyr::mutate(patientID=patient)
    }

    ctree_init = list("CCF_table"=CCF_table,
                      "drivers_table"=drivers_table,
                      "samples"=samples,
                      "patient"=patient)
    return(ctree_init)
}

if ( grepl(".rds\$", tolower("$ctree_input")) ) {
    best_fit = readRDS("$ctree_input")
    do_fit = TRUE

    if (!"data" %in% names(best_fit)) { best_fit[["data"]] =  data.frame(row.names=1:best_fit[["N"]])}

    ## viber
    if (class(best_fit) == "vb_bmm") {
        fn_name = VIBER::get_clone_trees
        subclonal_tool = "VIBER"

    # If only one sample inference ctree is not working with VIBER
    if (ncol(best_fit[["x"]]) <= 2) {
        cli::cli_alert_warning("Object of class {class(best_fit)} with only one sample is not supported by the function `get_clone_trees`!")
        do_fit = FALSE
    }

    best_fit[["data"]] = add_dummy_driver(best_fit[["data"]], variant_colname="gene", is_driver_colname="driver")
    }

    ## mobster
    if (class(best_fit) == "dbpmm") {
        fn_name = mobster::get_clone_trees
        subclonal_tool = "mobster"
        sample_id = unique(best_fit[["data"]][["sample_id"]])
        outdir = paste0(sample_id, "/")
        best_fit[["data"]] = add_dummy_driver(best_fit[["data"]], variant_colname="driver_label", is_driver_colname="is_driver")
    }

    if (class(best_fit) %in% c("vb_bmm", "dbpmm") & do_fit) {
        # VIBER or MOBSTER object
        trees = fn_name(x = best_fit)
    } else {
        cli::cli_alert_warning("Object of class {class(best_fit)} not supported.")
        do_fit = FALSE
    }

} else {
    do_fit = TRUE
    subclonal_tool = "pyclonevi"
    input_table = read.csv("$ctree_input", sep="\t")
    data_ctree = initialize_ctree_obj_pyclone(input_table)
    trees = ctrees(CCF_clusters = data_ctree[["CCF_table"]],
                   drivers = data_ctree[["drivers_table"]],
                   samples = data_ctree[["samples"]],
                   patient = data_ctree[["patient"]])
}

if (do_fit & !is.null(trees)) {
    ctree_output = paste0("ctree_", subclonal_tool)
    if (outdir != "" & !dir.exists(outdir)) dir.create(outdir, recursive=T)

    # plot the best tree
    top_phylo = plot(trees[[1]])

    # save rds and plots
    saveRDS(object=trees, file=paste0(outdir, opt[["prefix"]], "_", ctree_output, ".rds"))
    saveRDS(object=top_phylo, file=paste0(outdir, opt[["prefix"]], "_", ctree_output, "_plots.rds"))

    # save report plot
    phylos = ggplot2::ggplot()
    if (length(trees) > 1) phylos = lapply(trees[2:min(length(trees), 3)], plot)
    phylos = ggpubr::ggarrange(plotlist=phylos)

    ccf = ctree::plot_CCF_clusters(trees[[1]])
    info_transfer = ctree::plot_information_transfer(trees[[1]])
    clone_size = ctree::plot_clone_size(trees[[1]])

    report_fig = ggpubr::ggarrange(plotlist=list(ccf, info_transfer, top_phylo, clone_size, phylos), nrow=3, ncol=2)

    saveRDS(object=report_fig, file=paste0(outdir, opt[["prefix"]], "_", ctree_output, "_report.rds"))
    ggplot2::ggsave(plot=report_fig, filename=paste0(outdir, opt[["prefix"]], "_", ctree_output, "_report.pdf"), height=297, width=210, units="mm", dpi=200)
    ggplot2::ggsave(plot=report_fig, filename=paste0(outdir, opt[["prefix"]], "_", ctree_output, "_report.png"), height=297, width=210, units="mm", dpi=200)
}


# Version export

f = file("versions.yml","w")
ctree_version = sessionInfo()\$otherPkgs\$ctree\$Version
mobster_version = sessionInfo()\$otherPkgs\$mobster\$Version
viber_version = sessionInfo()\$otherPkgs\$VIBER\$Version
cli_version = sessionInfo()\$otherPkgs\$cli\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
tidyr_version = sessionInfo()\$otherPkgs\$tidyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
ggpubr_version = sessionInfo()\$otherPkgs\$ggpubr\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    ctree:", ctree_version), f)
writeLines(paste("    mobster:", mobster_version), f)
writeLines(paste("    VIBER:", viber_version), f)
writeLines(paste("    cli:", cli_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    tidyr:", tidyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    ggpubr:", ggpubr_version), f)
close(f)
