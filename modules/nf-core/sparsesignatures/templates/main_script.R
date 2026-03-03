#!/usr/bin/env Rscript

pkgs <- c("SparseSignatures", "dplyr", "tidyr", "tibble", "purrr", "stringr", "ggplot2", "patchwork")
sapply(pkgs, require, character.only = TRUE)

parse_args <- function(x) {
  # Remove brackets
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)

  # Split into key:value pairs
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))

  # Convert to tibble and split key:value
  tibble(arg = args_list) %>%
    separate(arg, into = c("key", "value"), sep = ":", extra = "merge", fill = "right") %>%
    dplyr::mutate(across(everything(), ~str_trim(.))) %>%
    dplyr::filter(!is.na(key) & !is.na(value)) %>%
    tibble::deframe()
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    genome = "${genome}",
    K = "2:10",
    nmf_runs = "10",
    iterations = "30",
    max_iterations_lasso = "10000",
    num_processes = "all",
    cross_validation_entries = "0.01",
    cross_validation_repetitions = "50",
    cross_validation_iterations = "5",
    lambda_values_alpha = "c(0.00, 0.01, 0.05, 0.10)",
    lambda_values_beta = "c(0.01, 0.05, 0.1, 0.2)",
    lambda_rate_alpha = "0",
    seed = "NULL"
)

args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]


# Script #

num_procs_string <- opt[["num_processes"]]

if (is.null(num_procs_string)) {
    stop("Missing required option: num_processes")
} else if (num_procs_string == "all") {
    n_procs <- Inf
} else {
    n_procs <- as.integer(opt[["num_processes"]])
}


patients_tsv = strsplit("$tsv_join", " ")[[1]]
tables = lapply(patients_tsv, FUN = function(p_table){
    read.delim(p_table, sep = "\\t", header=T) %>%
        mutate(across(everything(), as.character))
}
)
multisample_table = dplyr::bind_rows(tables)

#Extract input data information
input_data = multisample_table[,c("Indiv","chr","from","to","ref","alt")]
input_data = setNames(input_data, c("sample","chrom","start","end","ref","alt"))
input_data[["end"]] = input_data[["start"]]
input_data = input_data %>% mutate(start = as.integer(start), end = as.integer(end))

#Generate the patient vs mutation count matrix from mutation data
#Load reference human-genome specification.
#The user must select, among the available choices, the reference genome consistent with the mutation dataset.

load_genome = function(genome, input_data) {
    if (genome == "GRCh37") {
        library(BSgenome.Hsapiens.1000genomes.hs37d5)
        bsg = BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5

	if (any(grepl("^chr", input_data[["chrom"]]))) {
            input_data <- input_data %>% mutate(chrom = str_remove(chrom, "^chr"))
        }

    } else if (genome == "GRCh38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
        bsg = BSgenome.Hsapiens.UCSC.hg38

	if (all(!grepl("^chr", input_data[["chrom"]]))) {
            input_data <- input_data %>% mutate(chrom = paste0("chr", chrom))
        }

    }
    return(list(bsg = bsg, input_data = input_data))
}

data_list = load_genome(opt[["genome"]], input_data)
bsg = data_list[["bsg"]]
input_data = data_list[["input_data"]]

mut_counts = SparseSignatures::import.trinucleotides.counts(data=input_data, reference=bsg)

saveRDS(object = mut_counts, file = paste0(opt[["prefix"]], "_mut_counts.rds"))

# Load a reference SBS5 background signature from COSMIC
data(background)


seed_val <- if (!is.null(opt[["seed"]]) && opt[["seed"]] != "NULL") {
    as.integer(opt[["seed"]])
} else {
    42  # fallback default
}
set.seed(seed_val)

# Estimate the initial values of beta
starting_betas = SparseSignatures::startingBetaEstimation(x = mut_counts,
                                                        K = eval(parse(text=opt[["K"]])),
                                                        background_signature = background,
                                                        nmf_runs = as.integer(opt[["nmf_runs"]]),
                                                        seed = seed_val)

# Find the optimal number of signatures and sparsity level: rely on cross-validation
# higher number of CV repetitions corresponds to more accurate parameter estimates

cv_out = SparseSignatures::nmfLassoCV(
  x = mut_counts,
  K = eval(parse(text=opt[["K"]])),
  starting_beta = starting_betas,
  background_signature = background,
  nmf_runs = as.integer(opt[["nmf_runs"]]),
  lambda_values_alpha = eval(parse(text=opt[["lambda_values_alpha"]])),
  lambda_values_beta = eval(parse(text=opt[["lambda_values_beta"]])),
  cross_validation_entries = as.numeric(opt[["cross_validation_entries"]]),
  cross_validation_iterations = as.integer(opt[["cross_validation_iterations"]]),
  cross_validation_repetitions = as.integer(opt[["cross_validation_repetitions"]]),
  iterations = as.integer(opt[["iterations"]]),
  max_iterations_lasso = as.integer(opt[["max_iterations_lasso"]]),
  num_processes = n_procs,
  seed = seed_val
)

str(cv_out)


# Analyze the mean squared error results averaging over cross-validation repetitions
cv_mses = cv_out[["grid_search_mse"]][1,,]
cv_means_mse = matrix(sapply(cv_mses, FUN = mean),
                    nrow = dim(cv_mses)[1]
)

dimnames(cv_means_mse) = dimnames(cv_mses)

# Find the combination of parameters that yields the lowest MSE

print("CV MEAN MSE")
print(cv_means_mse)

print("MIN MEANS MSE")
print(min(cv_means_mse, na.rm = TRUE))

min_ii = which(cv_means_mse == min(cv_means_mse, na.rm = TRUE), arr.ind = TRUE)
min_Lambda_beta = rownames(cv_means_mse)[min_ii[1]]
min_Lambda_beta = as.numeric(gsub("_Lambda_Beta", "", min_Lambda_beta))
min_K = colnames(cv_means_mse)[min_ii[2]]
min_K = as.numeric(gsub("_Signatures", "", min_K))
best_params_config = data.frame(min_K, min_Lambda_beta)

saveRDS(object = cv_means_mse, file = paste0(opt[["prefix"]], "_cv_means_mse.rds"))
saveRDS(object = best_params_config, file = paste0(opt[["prefix"]], "_best_params_config.rds"))

# Discovering the signatures within the dataset: NMF Lasso
# Compute the signatures for the best configuration.

print(paste("MIN K =", min_K))

nmf_Lasso_out = SparseSignatures::nmfLasso(
  x = mut_counts,
  K = min_K,
  background_signature = background,
  lambda_rate_alpha = eval(parse(text=opt[["lambda_rate_alpha"]])),
  lambda_rate_beta = min_Lambda_beta,
  iterations = as.integer(opt[["iterations"]]),
  max_iterations_lasso = as.integer(opt[["max_iterations_lasso"]]),
  seed = seed_val
)

saveRDS(object = nmf_Lasso_out, file =  paste0(opt[["prefix"]], "_nmf_Lasso_out.rds"))

# Signature visualization
signatures = nmf_Lasso_out[["beta"]]
plot_signatures = SparseSignatures::signatures.plot(beta=signatures, xlabels=FALSE)

plot_exposure = nmf_Lasso_out[["alpha"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="PatientID") %>%
    tidyr::pivot_longer(cols=!"PatientID", names_to="Signatures", values_to="Exposures") %>%

    ggplot() +
    geom_bar(aes(x=PatientID, y=Exposures, fill=Signatures),
            position="stack", stat="identity") +
    theme(axis.text.x=element_text(angle=90,hjust=1),
            panel.background=element_blank(),
            axis.line=element_line(colour="black"))

plt_all = patchwork::wrap_plots(plot_exposure, plot_signatures, ncol=2) + patchwork::plot_annotation(title = "$meta.id")
ggplot2::ggsave(plot = plt_all, filename = paste0(opt[["prefix"]], "_plot_all.pdf"), width = 210, height = 297, units="mm", dpi = 200, device = "pdf")
saveRDS(object = plt_all, file = paste0(opt[["prefix"]], "_plot_all.rds"))

# version export
f = file("versions.yml","w")
SparseSignatures_version = sessionInfo()\$otherPkgs\$SparseSignatures\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
patchwork_version = sessionInfo()\$otherPkgs\$patchwork\$Version
bsg_hs37_version = sessionInfo()\$otherPkgs\$BSgenome.Hsapiens.1000genomes.hs37d5\$Version
bsg_hg38_version = sessionInfo()\$otherPkgs\$BSgenome.Hsapiens.UCSC.hg38\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    SparseSignatures:", SparseSignatures_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    patchwork:", patchwork_version), f)
writeLines(paste("    BSgenome.Hsapiens.1000genomes.hs37d5:", bsg_hs37_version), f)
writeLines(paste("    BSgenome.Hsapiens.UCSC.hg38", bsg_hg38_version), f)
close(f)
