#!/usr/bin/env Rscript

parse_args = function(x) {
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)
  # giving errors when we have lists like c(xxx, xxx) since it will separate it
  # args_list = unlist(strsplit(x, ', ')[[1]])
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
  # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
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

# Script
library(bascule)
library(ggplot2)
library(reticulate)
library(tidyverse)

untar("$counts_matrices", verbose=TRUE)
counts_folder = gsub(".tar.gz", "", "$counts_matrices")

counts_all = lapply(list.files(counts_folder, full.names=TRUE), function(file_id) {
  read.csv(file_id, sep="\t") %>% column_to_rownames(var="MutationType") %>%
    t() %>% as.data.frame()
}) %>% setNames(list.files(counts_folder))

n_types = length(counts_all)
types_names = names(counts_all)

# keep only samples present in all events types
samples_list = lapply(types_names, function(i) {
  tibble(counts_id=i, samples=rownames(counts_all[[i]]))
}) %>% bind_rows() %>%
  group_by(samples) %>%
  summarise(n_types=n()) %>%
  filter(n_types == !!n_types) %>%
  pull(samples) %>% as.character()

counts = lapply(counts_all, function(count_i) count_i[samples_list,])

reference_cat = eval(parse(text=opt[["reference_cat"]]))

# reorder columns
counts = lapply(types_names, function(type_id) {
  if (!is.null(reference_cat[[type_id]]))
    counts[[type_id]][,colnames(reference_cat[[type_id]])]
}) %>% setNames(types_names)


# set up reticulate
## define micromamba binary as main conda binary
options(reticulate.conda_binary=system("which micromamba", intern=TRUE))
## define the python binary - not necessary if RETICULATE_PYTHON is defined
# use_python(system("which python", intern=TRUE), required=TRUE)
## load pybascule package
py = import("pybascule")

x = fit(counts=counts,
        k_list=eval(parse(text=opt[["k_list"]])),
        reference_cat=reference_cat,
        keep_sigs=eval(parse(text=opt[["keep_sigs"]])),
        store_fits=TRUE,
        py=py)

x = refine_denovo_signatures(x)

clustering = FALSE
if (length(samples_list) > 1) {
  clustering = TRUE
  x = fit_clustering(x, cluster=as.integer(opt[["cluster"]]))
  x = merge_clusters(x, cutoff=as.numeric(opt[["cutoff"]]))
}

x = convert_dn_names(x, cutoff=as.numeric(opt[["cutoff"]]))

# plots
pl_signatures = plot_signatures(x)
pl_exposures = plot_exposures(x)
pl_centroids = plot_centroids(x)

if (clustering) {
  figure = patchwork::wrap_plots(pl_signatures,
                                 pl_exposures + theme(legend.position="bottom"),
                                 pl_centroids + theme(legend.position="bottom"),
                                 design="aaaa\nbbbc")
} else {
  figure = patchwork::wrap_plots(pl_signatures,
                                 pl_exposures + theme(legend.position="bottom"),
                                 ncol=1)
}

# save objects
saveRDS(x, file=paste0(opt[["prefix"]], "_fit_bascule.rds"))
saveRDS(figure, paste0(opt[["prefix"]], "_plots_bascule.rds"))
ggplot2::ggsave(plot=figure, filename=paste0(opt[["prefix"]], "_plots_bascule.pdf"), height=210, width=210, units="mm")
ggplot2::ggsave(plot=figure, filename=paste0(opt[["prefix"]], "_plots_bascule.png"), height=210, width=210, units="mm", dpi=200)


# version export
f = file("versions.yml","w")
bascule_version = sessionInfo()\$otherPkgs\$bascule\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
reticulate_version = sessionInfo()\$otherPkgs\$reticulate\$Version
tidyverse_version = sessionInfo()\$otherPkgs\$tidyverse\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    bascule:", bascule_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    reticulate:", reticulate_version), f)
writeLines(paste("    tidyverse:", tidyverse_version), f)
close(f)
