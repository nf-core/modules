#!/usr/bin/env Rscript

# Differential RNA-editing test across sample groups, using the beta-binomial
# (or binomial / quasibinomial / betareg-like) GLM from Bullseye:
# https://github.com/mflamand/Bullseye (example_data/markdown/functions.R, MIT License,
# Copyright (c) 2021 mflamand). The upstream functions are reproduced below until
# Bullseye ships them as a package.

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(aod)
})

parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })
    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[!is.na(parsed_args)]
}

opt <- list(
    prefix         = ifelse('${task.ext.prefix}' == 'null', '${meta.id}', '${task.ext.prefix}'),
    coverage       = '${coverage}',
    mutations      = '${mutations}',
    samplesheet    = '${samplesheet}',
    sample_column  = 'sample',
    reference_level = '',
    design         = '~condition',
    min_cov        = 10,
    min_edit       = 0.05,
    link           = 'betabin',
    cores          = as.integer('${task.cpus}')
)
opt_types <- lapply(opt, class)

args_opt <- parse_args('${task.ext.args}')
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    }
    args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    opt[[ao]] <- args_opt[[ao]]
}

#######################################################################
## Bullseye GLM functions (upstream, unchanged except where noted)   ##
#######################################################################

### Function for betabinomial distribution
glm_comp_bbn<-function(df, design, lev,comp){

  facts_comp<-levels(df[[comp]])[[1]]
  df<-subset(df, df\$cov != 0)
  df\$non_mut<-df\$cov-df\$mut

  # to avoid issues when no mutations are found in any samples, we add pseudocounts to the model. These are 0.5% of the geometric mean of the coverage data.
  baseMean<-(df %>% summarise(baseMean= exp(mean(log(cov+0.5)))*0.005))\$baseMean
  df<-df%>% mutate(cov= cov +baseMean,mut=mut+baseMean)

  facts<-levels(droplevels(df[[comp]]))
  test_res<-NULL

  design1<-design[[2]]
  lev=paste0(design1,lev)
  levLFC<-paste0(lev,"LFC")

  if (length(facts) < 2  || !facts_comp %in% facts){
    p<-as.vector(rep(NA,length(lev)),mode="list")
    names(p)<-lev

  }
  else{
    suppressWarnings(mod<-aod::betabin(formula(paste0("cbind(mut,cov-mut)",deparse(design))), random = ~1, data=df))
    test_res<-try(vcov(mod),silent=TRUE)
  }

  # Patched: upstream has a stray bare "." statement here that always throws
  # "object '.' not found".
  if(is.null(test_res) ){
    p<-as.vector(rep(NA,length(lev)),mode="list")
    names(p)<-lev
    return(p)
  }
  else{
    p<-vector("list", length(lev))
    n<-vector("list", length(lev))
    LFC<-as.vector(rep(0,length(lev)),mode="list")
    names(LFC)<-levLFC

    okbeta <- !is.na(coef(mod, na.rm = FALSE))
    aa <- coef(mod)[okbeta]

    index <- which(names(aa) %in% lev)

    for(x in index){

      b<-coef(mod)
      V<-vcov(mod)
      L <- matrix(rep(0, length(b)), ncol = length(b))
      L[x] <- 1
      H0 <- rep(0, nrow(L))
      f <- L %*% b

      mat <- qr.solve(L %*% V %*% t(L))
      stat <- t(f - H0) %*% mat %*% (f - H0)
      # Patched: upstream's 1 - pchisq() underflows to exactly 0 for large statistics.
      p[x] <- pchisq(stat, df = nrow(L), lower.tail = FALSE)
      n[x] <- names(b)[x]
    }

    names(p)<- n

    for ( x in lev) {
      if (!x %in% names(p)){
        p[[x]]<- NA   # if it was missing because of a lack of coverage, return NA as p-value
      }
    }
  }
  return(p[ names(p) %in% lev])
}

# Function for GLM binomial, quasibinomial or betargression (like) distribution
glm_comp_bn<-function(df, design, lev,comp,link="binomial"){

  facts_comp<-levels(df[[comp]])[[1]]
  df<-subset(df, df\$cov != 0)

  baseMean<-(df %>% summarise(baseMean= exp(mean(log(cov+0.5)))*0.005))\$baseMean
  df<-df%>% mutate(cov= cov +baseMean,mut=mut+baseMean)

  canonicalOrder <- function(term) {
    tt <- strsplit(term, ":")
    tt <- lapply(tt, sort)
    sapply(tt, paste, collapse = ":")
  }

  explicit1<-function (formula) {
    if (length(formula) == 1)
      return(formula == 1)
    if (!(formula[[1]] == "+" || formula[[1]] == "*" || formula[[1]] ==
          "/" || formula[[1]] == "^" || formula[[1]] == "~"))
      return(FALSE)
    if (length(formula) == 3) {
      (formula[[2]] == 1) || explicit1(formula[[2]]) || explicit1(formula[[3]])
    }
    else {
      (formula[[2]] == 1) || explicit1(formula[[2]])
    }
  }

  facts<-levels(droplevels(df[[comp]]))
  if (length(facts) < 2  || !facts_comp %in% facts){
    p<-as.vector(rep(1,length(lev)),mode="list")
    LFC<-as.vector(rep(0,length(lev)),mode="list")
    design<-design[[3]]
    lev=paste0(design,lev)
    levLFC<-paste0(lev,"LFC")
    names(p)<-lev
    names(LFC)<-levLFC
    lev<-c(lev,levLFC)
    p<-c(p,LFC)
    return(p[ names(p) %in% lev])
  }
  else{
    mod<-glm(formula=design,data=df,weights = cov, family=link)
    design<-design[[3]]
    lev=paste0(design,lev)
    levLFC<-paste0(lev,"LFC")

    p<-vector("list", length(lev))
    LFC<-as.vector(rep(0,length(lev)),mode="list")
    names(LFC)<-levLFC
    n<-vector("list", length(lev))

    design<-formula(paste0("~",design))

    if (inherits(design, "formula")) {
      test_intercept <- explicit1(design)
      design <- attr(terms(design), "term.labels")
    }
    else test_intercept <- FALSE

    okbeta <- !is.na(coef(mod, na.rm = FALSE))
    tt <- attr(terms(mod), "term.labels")
    aa <- attr(model.matrix(mod), "assign")[okbeta]

    index <- which(aa %in% match(canonicalOrder(design), canonicalOrder(tt)))

    if (test_intercept) {
      if (attr(terms(model), "intercept"))
        index <- unique(c(1, index))
      else stop("model does not have an intercept")
    }

    for(x in index){
      beta <- coef(mod)[x]
      name_fact<-(names(beta))
      V<-NULL
      V<-try(vcov(mod)[x, x],silent=TRUE)
      if(is.null(V) || is.na(V) || class(V)=="try-error" || V == 0){
        p[x]<-1
      }
      else{
        chisq <- beta %*% solve(V) %*% beta
        pval = pchisq(chisq, length(2), lower.tail = FALSE)
        p[x]<-pval[1,1] # get p-value from result
      }
      n[x]<-name_fact
      LFC[x]<-beta
    }
    names(p)<- n

    for ( x in lev) {
      if (!x %in% names(p)){
        p[[x]]<- NA   # if it was missing because of a lack of coverage, return NA as p-value
        LFC[[paste0(x,"LFC")]]<- 0
      }
    }
  }
  p<-c(p,LFC)
  lev<-c(lev,levLFC)
  return(p[ names(p) %in% lev])

}

#function to calculate glm and p value
Bullseye<-function(data_cov,data_mut, colData=NULL, filter.col=NULL, filter.with=NULL, min=0.05, filter.sites=NULL,dataset=NULL,design=NULL, min.cov=10, return.df.only=FALSE, link='betabin',numCores=1){

  #filter.sites : needs to be "all" or "called", if all: will only keep sites with at least min edit in one or the other samples. If "called": will keep samples with at least min in both conditions.
  if(link != 'betabin' && link != 'bin' && link != 'quasibinomial' && link !="betareg") stop('Please use link="betabin", link="bin", link="quasibinomial" or link="betareg"')

  if(is.null(colData)) stop("please provide colData with sample information for filtering")
  if(is.null(design)) stop("please provide design formula with a factor found in colData (eg. ~Genotype)")

  if(!is.null(filter.col)){
   colData<- colData[colData[[filter.col]]==filter.with,] #base version
  }

  if(!is.factor(colData[[ design[[2]] ]])){
    colData[[ design[[2]] ]]<- as.factor(colData[[ design[[2]] ]])
  }

  colData[[ design[[2]] ]]<- droplevels(colData[[ design[[2]] ]])

  design_levels<-levels(colData[[ design[[2]] ]])

  #only keep matrix columns from columns in colData

  if(! "n" %in% colnames(data_cov))
  {
    data_cov\$n = 1
    data_mut\$n = 1
  }

  data_cov<-data_cov[,c("n",rownames(colData))]
  data_mut<-data_mut[,c("n",rownames(colData))]

  # We can filter based on a list of sites found in .dataset. Else we can filter and keep a subset of samples from colData.
  if(!is.null(dataset)){
    dataset\$site<- paste(dataset\$chr,dataset\$start,dataset\$end,dataset\$strand,sep="_")
    data_cov<-data_cov[rownames(data_cov) %in% dataset\$site,]
    data_mut<-data_mut[rownames(data_mut) %in% dataset\$site,]
  }

  n_site<-data_cov\$n

  # filter based on min coverage, replace with 0 to exclude from analysis
  data_cov[data_cov < min.cov*data_cov\$n] <- 0
  data_mut[data_cov < min.cov*data_mut\$n] <- 0

  data_cov\$n<-n_site
  data_mut\$n<-n_site
  # calculate editing rate in each sample
  data_edit<-(data_mut/data_cov)*n_site
  data_edit[is.na(data_edit)]<-NA

  # first we get the average editing in called sites
  data_edit_perc<-data_edit
  data_edit_perc[data_edit_perc< min]<-NA
  data_edit_perc[is.na(data_edit_perc)]<-NA  # NaN to NA if any

  data_cov <- data_cov %>% rownames_to_column("site") # change to tibble style for pivoting and nesting of dataframe
  data_mut <- data_mut %>% rownames_to_column("site")
  data_edit <- data_edit %>% rownames_to_column("site")
  data_edit_perc <- data_edit_perc %>% rownames_to_column("site")

  L_cov<- data_cov %>%
    pivot_longer(names_to = "sample",cols = !c(site,n) ,values_to = "cov") %>%
    left_join(., colData %>% rownames_to_column("sample"), by="sample")

  L_mut<-data_mut %>%
    pivot_longer(names_to = "sample",cols = !c(site,n),values_to = "mut")

  data_long <- suppressMessages(dplyr::left_join(L_cov,L_mut))
  data_long<-data_long %>% mutate(edit=mut/cov)
  n_data<-data_long%>% group_by(site) %>% nest()

  if (link == "bin"){
    results<-parallel::mclapply(X = n_data\$data, FUN = glm_comp_bn,design=formula(paste0("mut/cov",deparse(design))), lev=design_levels,comp=design[[2]],link="binomial",mc.cores = numCores)
  }
  else if (link =="quasibinomial"){
    results<-parallel::mclapply(X = n_data\$data, FUN = glm_comp_bn,design=formula(paste0("mut/cov",deparse(design))), lev=design_levels,comp=design[[2]],link="quasibinomial",mc.cores = numCores)
  }
  else if (link =="betabin"){
    results<-parallel::mclapply(X = n_data\$data, FUN = glm_comp_bbn,design=design, lev=design_levels,comp=design[[2]],mc.cores = numCores)
  }
  else if (link =="betareg"){
    results<-parallel::mclapply(X = n_data\$data, FUN = glm_comp_bn,design=formula(paste0("mut/cov",deparse(design))), lev=design_levels,comp=design[[2]],link=quasi(link="logit", variance="mu(1-mu)"),mc.cores = numCores)
  }
  results<-do.call(rbind.data.frame, results)
  site<- n_data\$site
  results<-cbind(site,results)

  # first level is always first and is fixed.
  first.level<-design_levels[[1]]
  first.level.regex<-paste0(".*",first.level,".*")

  data_edit_A<-data.frame("site"='')

  i<-length(design_levels)
  while( i >= 2){

    second.level<-design_levels[[i]]
    second.level.regex<-paste0(".*",second.level,".*")

    data_edit_A<-suppressMessages(dplyr::mutate(.data=data_cov,
                                                !!paste0(first.level,"_sum_cov") := rowSums(select(data_cov, matches( first.level.regex,perl=TRUE)),na.rm=TRUE)) %>%
                                    dplyr::mutate(!!paste0(first.level,"_mut") := rowSums(select(data_mut, matches(first.level.regex,perl=TRUE)),na.rm=TRUE)) %>%
                                    mutate(!!first.level := (.data[[paste0(first.level,"_mut")]]/.data[[paste0(first.level,"_sum_cov")]]) *n_site) %>%
                                    mutate(!!paste0(first.level,"_min_cov") := rowSums(!is.na(select(.data=data_edit, matches(first.level.regex,perl=TRUE)))),
                                           !!paste0(first.level,"_called_sites") := rowSums((select(.data=data_edit, matches(first.level.regex,perl=TRUE)))   >= min,na.rm=TRUE )) %>%
                                    mutate(!!paste0(first.level,"_average_called") := rowMeans(select(.data=data_edit_perc, matches(first.level.regex,perl=TRUE)),na.rm = TRUE)) %>%
                                    dplyr::mutate(!!paste0(second.level,"_sum_cov") := rowSums(select(.data=data_cov, matches( second.level.regex,perl=TRUE)),na.rm = TRUE)) %>%
                                    dplyr::mutate(!!paste0(second.level,"_mut") := rowSums(select(.data=data_mut, matches(second.level.regex,perl=TRUE)),na.rm = TRUE)) %>%
                                    # Note (not patched): unlike glm_comp_bbn()'s p-value, which adds a
                                    # pseudocount before fitting, this descriptive log2FoldChange is computed
                                    # directly from the raw pooled editing rates. Whenever the first (reference)
                                    # level has zero mutations at a site, this divides by zero and reports Inf
                                    # rather than a finite fold-change - seen on ~26% of sites in the module's
                                    # own test data. The p-value/padj columns are unaffected and remain valid.
                                    mutate(!!second.level := (.data[[paste0(second.level,"_mut")]]/.data[[paste0(second.level,"_sum_cov")]] *n_site),
                                           !!paste0(second.level,"_log2FoldChange") := log2(.data[[second.level]]/.data[[first.level]])) %>%
                                    mutate(!!paste0(second.level,"_min_cov") := rowSums(!is.na(select(.data=data_edit, matches(second.level.regex,perl=TRUE)))),
                                           !!paste0(second.level,"_called_sites") := rowSums((select(.data=data_edit, matches(second.level.regex,perl=TRUE))) >= min,na.rm=TRUE)) %>%
                                    mutate(!!paste0(second.level,"_average_called") := rowMeans(select(.data=data_edit_perc, matches(second.level.regex,perl=TRUE)),na.rm = TRUE)) %>%
                                    left_join(.,data_edit_A) %>% left_join(.,results %>% select(site, !!paste0(second.level,"_pvalue"):=matches( paste0(".*",second.level,"\$")),
                                                                                                !!paste0(second.level,"_LFC"):=matches(paste0(".*",second.level,"LFC") )) ) )

    i<-i-1
  }

  data_edit_A<- data_edit_A %>% mutate(baseEdit = rowSums(select(., matches("average_called")),na.rm = TRUE)*100*n/rowSums(select(., matches("called_sites")),na.rm = TRUE)) %>% mutate(baseEdit=replace_na(baseEdit,0))
  data_edit_A<-data_edit_A %>% dplyr::select(-rownames(colData))

  if (!is.null(filter.sites)){
    # Patched: upstream uses R.utils::printf, an otherwise unneeded dependency.
    message(sprintf("Filtering sites with less than %.1f %% editing", min*100))
    if (filter.sites == "called"){
      data_edit_A<-data_edit_A %>%
        dplyr::filter(baseEdit >= min)
    }
    else if (filter.sites == "all"){
      data_edit_A<-data_edit_A %>%
        dplyr::filter(.data[[first.level]] >= min | .data[[second.level]] >= min)
    }
  }

  return(data_edit_A)

}

#######################################################################
## Inputs                                                            ##
#######################################################################

cov <- fread(opt\$coverage)
mut <- fread(opt\$mutations)
observations <- as.data.frame(fread(opt\$samplesheet))
design <- as.formula(opt\$design)
group_col <- all.vars(design)[1]

if (!opt\$sample_column %in% colnames(observations)) {
    stop(sprintf("samplesheet must contain the sample column '%s'", opt\$sample_column))
}
if (!group_col %in% colnames(observations)) {
    stop(sprintf("samplesheet must contain the design column '%s'", group_col))
}

site_cols <- c("#chr", "start", "end", "strand")
if (!all(site_cols %in% colnames(cov)) || !identical(colnames(cov), colnames(mut))) {
    stop("coverage and mutation matrices must have identical columns, including ",
         paste(site_cols, collapse = ", "), ", plus one column per sample")
}
site_id <- paste(cov[["#chr"]], cov\$start, cov\$end, cov\$strand, sep = "_")
if (!identical(site_id, paste(mut[["#chr"]], mut\$start, mut\$end, mut\$strand, sep = "_"))) {
    stop("coverage and mutation matrices must list the same sites in the same order")
}

samples <- as.character(observations[[opt\$sample_column]])
missing_samples <- setdiff(samples, colnames(cov))
if (length(missing_samples) > 0) {
    stop("samplesheet samples not found in the matrix columns: ", paste(missing_samples, collapse = ", "))
}

groups <- as.character(observations[[group_col]])
if (length(unique(groups)) < 2) {
    stop(sprintf("design column '%s' needs at least 2 groups, found: %s",
                 group_col, paste(unique(groups), collapse = ", ")))
}
if (opt\$reference_level != '' && !opt\$reference_level %in% groups) {
    stop(sprintf("reference_level '%s' not found in design column '%s' (levels: %s)",
                 opt\$reference_level, group_col, paste(unique(groups), collapse = ", ")))
}

# Bullseye() assigns sample columns to groups by regex (".*<group>.*"), so name
# each column <group>__<sample> and reject labels that also match another group.
new_names <- paste0(groups, "__", samples)
for (lev in unique(groups)) {
    hit <- grepl(paste0(".*", lev, ".*"), new_names, perl = TRUE)
    if (!identical(hit, groups == lev)) {
        stop(sprintf("group '%s' also matches samples of another group (%s); group labels must not occur inside other group or sample names",
                     lev, paste(new_names[hit & groups != lev], collapse = ", ")))
    }
}

#######################################################################
## Run                                                               ##
#######################################################################

outfile <- paste0(opt\$prefix, ".glm.tsv")

if (nrow(cov) == 0) {
    warning("coverage/mutation matrices contain no sites; writing an empty result")
    fwrite(data.table(site = character()), file = outfile, sep = "\t")
} else {
    # rownames<- is dropped by data.table, hence data.frame; Bullseye() re-adds
    # rownames as a 'sample' column, so colData must not carry one.
    data_cov <- as.data.frame(cov)[, samples, drop = FALSE]
    data_mut <- as.data.frame(mut)[, samples, drop = FALSE]
    colnames(data_cov) <- new_names
    colnames(data_mut) <- new_names
    rownames(data_cov) <- site_id
    rownames(data_mut) <- site_id
    col_data <- observations[, setdiff(colnames(observations), opt\$sample_column), drop = FALSE]
    rownames(col_data) <- new_names
    # Bullseye() compares every level against the first one
    if (opt\$reference_level != '') {
        col_data[[group_col]] <- relevel(factor(col_data[[group_col]]), ref = opt\$reference_level)
    }

    result <- Bullseye(
        data_cov       = data_cov,
        data_mut       = data_mut,
        colData        = col_data,
        design         = design,
        min            = opt\$min_edit,
        min.cov        = opt\$min_cov,
        return.df.only = TRUE,
        link           = opt\$link,
        numCores       = opt\$cores
    )
    result <- as.data.frame(result)
    for (pcol in grep("_pvalue\$", colnames(result), value = TRUE)) {
        padj <- sub("_pvalue\$", "_padj", pcol)
        result[[padj]] <- p.adjust(result[[pcol]], method = "BH")
        pos <- match(pcol, colnames(result))
        result <- result[, c(colnames(result)[seq_len(pos)], padj, setdiff(colnames(result)[-seq_len(pos)], padj))]
    }
    fwrite(as.data.table(result), file = outfile, sep = "\t")
}

sink(paste0(opt\$prefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', as.character(getRversion())),
        paste('    r-aod:', as.character(packageVersion('aod')))
    ),
    'versions.yml'
)
