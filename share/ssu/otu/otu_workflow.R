load(here("page_build", "otu_part1.rdata"))

set.dir(output=pipelineFiles/)

make.file(inputdir=$DATA, type=$TYPE, prefix=shrimp)
# OR
make.file(inputdir=01_TRIMMED_DATA/, type=gz, prefix=shrimp)

make.contigs(file=current, processors=$PROC)
# OR
make.contigs(file=shrimp.files, processors=30)

summary.seqs(fasta=current, count=current, processors=$PROC)
summary.seqs(fasta=shrimp.trim.contigs.fasta, count=shrimp.contigs.count_table, processors=30)

count.groups(count=current)
count.groups(count=shrimp.contigs.count_table)

screen.seqs(fasta=current, count=current, maxambig=0, minlength=252, maxlength=254, maxhomop=6, processors=$PROC)
screen.seqs(fasta=current, count=current, maxambig=0, minlength=252, maxlength=254, maxhomop=6, processors=30)

summary.seqs(fasta=current, count=current, processors=$PROC)
summary.seqs(fasta=shrimp.trim.contigs.good.fasta, count=shrimp.contigs.good.count_table, processors=30)

count.groups(count=current)
count.groups(count=shrimp.contigs.good.count_table)

unique.seqs(fasta=current, count=current)
unique.seqs(fasta=shrimp.trim.contigs.good.fasta, count=shrimp.contigs.good.count_table)

summary.seqs(count=current, processors=$PROC)
summary.seqs(count=shrimp.trim.contigs.good.count_table, processors=30)

pcr.seqs(fasta=$REF_LOC/silva.seed_v138_2.align, start=13862, end=23445, keepdots=F, processors=$PROC)
pcr.seqs(fasta=reference_dbs/silva.seed_v138_2.align, start=13862, end=23445, keepdots=F, processors=30)

rename.file(input=pipelineFiles/silva.nr_v132.pcr.align, new=pipelineFiles/$ALIGNREF)
rename.file(input=pipelineFiles/silva.seed_v138_2.pcr.align, new=pipelineFiles/silva.v4.fasta)

summary.seqs(fasta=pipelineFiles/$ALIGNREF, processors=$PROC)
summary.seqs(fasta=pipelineFiles/silva.v4.fasta, processors=30)

align.seqs(fasta=current, reference=$REF_LOC/$ALIGNREF, processors=$PROC)
align.seqs(fasta=shrimp.trim.contigs.good.unique.fasta, reference=pipelineFiles/silva.v4.fasta, processors=30)

summary.seqs(fasta=current, count=current, processors=$PROC)
summary.seqs(fasta=shrimp.trim.contigs.good.unique.align, count=shrimp.trim.contigs.good.count_table, processors=30)

screen.seqs(fasta=current, count=current, start=1968, end=11550, processors=$PROC)
screen.seqs(fasta=shrimp.trim.contigs.good.unique.align, count=shrimp.trim.contigs.good.count_table, start=1, end=9583, processors=30)

summary.seqs(fasta=current, count=current, processors=$PROC)
summary.seqs(fasta=shrimp.trim.contigs.good.unique.good.align, count=shrimp.trim.contigs.good.good.count_table, processors=30)

count.groups(count=current)
count.groups(count=shrimp.trim.contigs.good.good.count_table)

filter.seqs(fasta=current, vertical=T, trump=., processors=$PROC)
filter.seqs(fasta=shrimp.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=30)

unique.seqs(fasta=current, count=current)
unique.seqs(fasta=shrimp.trim.contigs.good.unique.good.filter.fasta, count=shrimp.trim.contigs.good.good.count_table)

summary.seqs(fasta=current, count=current, processors=$PROC)
summary.seqs(fasta=shrimp.trim.contigs.good.unique.good.filter.unique.fasta, count=shrimp.trim.contigs.good.unique.good.filter.count_table, processors=30)

system(cp pipelineFiles/shrimp.trim.contigs.good.unique.good.filter.unique.fasta pipelineFiles_med/)
system(cp pipelineFiles/shrimp.trim.contigs.good.unique.good.filter.count_table pipelineFiles_med/)

pre.cluster(fasta=current, count=current, diffs=2, processors=$PROC)

summary.seqs(fasta=current, count=current, processors=$PROC)

count.groups(count=current)

# tmp_accnos <- readr::read_delim(here(work_here, "nc_screen/shrimp.files"),
#                                 delim = "\t", col_names = FALSE)
# tmp_accnos[, 2:3] <- NULL
# tmp_accnos <- tmp_accnos[grepl("Control_", tmp_accnos$X1), ]
# readr::write_delim(tmp_accnos, file = here(work_here, "nc_screen/nc_samples.accnos"),
#                    col_names = FALSE)

get.groups(fasta=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.count_table, accnos=nc_samples.accnos)

rename.file(input=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, new=nc.fasta)
rename.file(input=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, new=nc.count_table)

summary.seqs(fasta=nc.fasta, count=nc.count_table, processors=$PROC)

list.seqs(count=nc.count_table)

count.seqs(count=nc.count_table, compress=f)

get.seqs(accnos=nc.accnos, fasta=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

rename.file(input=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, new=subset.fasta)
rename.file(input=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, new=subset.count_table)

count.seqs(count=subset.count_table, compress=f)

# full_count_tab <- readr::read_delim(here(work_here, "nc_screen/subset.full.count_table"),
#                                     delim = "\t", col_names = TRUE)
# # figure out which columns to use
# control_cols     <- grep("^Control_", names(full_count_tab), value = TRUE)
# noncontrol_cols  <- setdiff(names(full_count_tab)[-(1:2)], control_cols)
# 
# # now do the rowwise sums
# read_totals <- full_count_tab %>%
#   rowwise() %>%
#   mutate(
#     total_reads_nc   = sum(c_across(all_of(control_cols)), na.rm = TRUE),
#     total_reads_non_nc = sum(c_across(all_of(noncontrol_cols)), na.rm = TRUE)
#   ) %>%
#   ungroup() %>%
#   select(1, 2, total_reads_nc, total_reads_non_nc)
# 
# read_totals <- read_totals %>% dplyr::rename("total_reads" = 2)

#head_read_totals <- head(read_totals)
#nc_dim <- dim(read_totals)[1]
pandoc.table(read_totals[1:6, 1:4], 
             emphasize.rownames = FALSE, 
             split.tables = Inf)

# tmp_read_totals <- read_totals %>%
#   dplyr::mutate(perc_reads_in_nc = 100*(
#     total_reads_nc / (total_reads_nc + total_reads_non_nc)),
#                 .after = "total_reads_non_nc")
# tmp_read_totals$perc_reads_in_nc <-
#   round(tmp_read_totals$perc_reads_in_nc, digits = 6)

# control_cols     <- grep("^Control_", names(full_count_tab), value = TRUE)
# noncontrol_cols  <- setdiff(names(full_count_tab)[-(1:2)], control_cols)
# # rowwise tally of non-zero columns
# samp_totals <- full_count_tab %>%
#   rowwise() %>%
#   mutate(
#     num_nc_samp     = sum(c_across(all_of(control_cols)) != 0, na.rm = TRUE),
#     num_non_nc_samp = sum(c_across(all_of(noncontrol_cols)) != 0, na.rm = TRUE)
#   ) %>%
#   ungroup() %>%
#   select(1, num_nc_samp, num_non_nc_samp)

# samp_totals$total_samp <- samp_totals$num_nc_samp + samp_totals$num_non_nc_samp
# samp_totals <- samp_totals %>%
#   dplyr::relocate("total_samp", .after = "Representative_Sequence")
# samp_totals <- samp_totals %>%
#   dplyr::mutate(perc_nc_samp =
#                   100*( num_nc_samp / (num_nc_samp + num_non_nc_samp)),
#                   .after = "num_non_nc_samp")

# nc_check <- dplyr::left_join(tmp_read_totals, samp_totals, by = "Representative_Sequence")
# write_delim(nc_check, here(work_here, "nc_screen/reads_in_nc_samples.txt"),
#     delim = "\t")

pandoc.table(nc_check[1:3, 1:9], 
             emphasize.rownames = FALSE, 
             split.tables = Inf)

# nc_remove <- nc_check %>%
#   filter(perc_reads_in_nc > 10 | perc_nc_samp > 10)

# nc_remain <- dplyr::anti_join(nc_check, nc_remove)
# 
# rem_nc_reads <- sum(nc_remove$total_reads_nc)
# rem_sam_reads <- sum(nc_remove$total_reads_non_nc)
# per_reads_rem <- round(100*( rem_nc_reads / (rem_nc_reads + rem_sam_reads)),
#                        digits = 3)
# 
# ret_nc_reads <- sum(nc_remain$total_reads_nc)
# ret_sam_reads <- sum(nc_remain$total_reads_non_nc)
# per_reads_ret <- round(100*( ret_nc_reads / (ret_nc_reads + ret_sam_reads)),
#                        digits = 3)

# write_delim(
#   data.frame(nc_remove$Representative_Sequence),
#   here(work_here, "nc_screen/nc_repseq_remove.accnos"),
#   col_names = FALSE)

remove.seqs(accnos=nc_repseq_remove.accnos, fasta=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.count_table)

count.groups(count=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table)

# tmp_before <- read_tsv(
#   here(work_here,
#        "nc_screen/shrimp.trim.contigs.good.unique.good.filter.unique.precluster.count.summary"),
#   col_names = FALSE,
#   col_select = 1
# )
# 
# tmp_after <- read_tsv(
#   here(work_here,
#        "nc_screen/shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.count.summary"),
#   col_names = FALSE,
#   col_select = 1
# )
# tmp_nc_lost <- anti_join(tmp_before, tmp_after)
# tmp_nc_lost$X1

# nc_to_remove <- semi_join(tmp_before, tmp_after)
# nc_to_remove <- nc_to_remove %>%
#   dplyr::filter(
#     stringr::str_starts(X1, "Control")
#     )
# readr::write_delim(nc_to_remove,
#                    file = here(work_here, "nc_screen/nc_samples_remove.accnos"),
#                    col_names = FALSE)

# nc_to_remove

remove.groups(count=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, fasta=shrimp.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, accnos=nc_samples_remove.accnos)

summary.seqs(fasta=current, count=current, processors=30)

count.groups(count=current)

chimera.vsearch(fasta=current, count=current, dereplicate=t, processors=$PROC)

summary.seqs(fasta=current, count=current, processors=30)

count.groups(count=current)

wget https://manichanh.vhir.org/gsrdb/GSR-DB_V4_cluster-1.tar.gz
tar -xvzf GSR-DB_V4_cluster-1.tar.gz

cp GSR-DB_V4_cluster-1_taxa.txt tmp0.txt
sed '1d' tmp0.txt > tmp1.txt

sed -E 's/s__.*//g' tmp1.txt > tmp2.txt
sed -E 's/[a-zA-Z]__//g' tmp2.txt > gsrdb.tax
cp GSR-DB_V4_cluster-1_seqs.fasta gsrdb.fasta

classify.seqs(fasta=current, count=current, reference=reference_dbs/gsrdb.fasta, taxonomy=reference_dbs/gsrdb.tax, processors=3)

remove.lineage(fasta=current, count=current, taxonomy=current, taxon=$CONTAMINENTS)

summary.tax(taxonomy=current, count=current)

count.groups(count=current)

# # this is to generate a read change table DO NOT RUN HERE
# # gather all count.summary file together in one directory
# library(dplyr)
# library(purrr)
# library(readr)
# files <- file.info(list.files("read_change_summary", pattern = "\\.summary$", full.names = TRUE))
# files <- files[order(files$mtime), ]
# files <- as.character(row.names(files))
# col_names <- c("SampleID",
#                "raw_reads",
#                "initial_screen",
#                "screen_after_align",
#                "remove_nc_reads",
#                "remove_nc",
#                "nochim",
#                "no_contam")
# 
# df <- files %>%
#   map(~ read_tsv(.x, col_names = c("SampleID", "value"), col_types = "ci")) %>%
#   reduce(left_join, by = "SampleID")
# 
# # Rename columns
# names(df) <- col_names
# 
# # Arrange rows (optional: if you want custom order of samples)
# # For example, put controls first:
# df <- df %>%
#   arrange(factor(SampleID, levels = sort(unique(SampleID))))
# # Write to tab-delimited file
# write_tsv(df, "mothur_pipeline_read_changes.txt")

# read_change <- read_tsv(
#   here(work_here, "mothur_pipeline_read_changes.txt"),
#   col_names = TRUE
# )

rename.file(fasta=current, count=current, taxonomy=current, prefix=final)

cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=4, cluster=f, processors=$PROC)
cluster.split(file=final.file, count=final.count_table, processors=$PROC)

system(mkdir cluster.split.gsrdb) 
system(mv final.opti_mcc.list cluster.split.gsrdb/) 
system(mv final.file cluster.split.gsrdb/) 
system(mv final.dist cluster.split.gsrdb/)

dist.seqs(fasta=final.fasta, cutoff=0.03, processors=\$PROC) cluster(column=final.dist, count=final.count_table)

# readr::write_delim(read_change,
#   here(share_here, "mothur_pipeline_read_changes.txt"),
#   col_names = TRUE, delim = "\t"
# )
# file.copy(here(work_here,  "nc_screen/reads_in_nc_samples.txt"),
#           here(share_here),
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)

# objects()
# gdata::keep(nc_check, nc_remain, nc_remove, nc_to_remove, nc_dim,
#             no_read_sum, per_reads_rem, per_reads_ret,
#             rem_nc_reads, rem_sam_reads, ret_nc_reads,
#             ret_sam_reads, read_change, nc_tmp1, nc_tmp2,
#             read_totals,
#             sure = TRUE)
# save.image(here("page_build", "otu_part1.rdata"))


# -------------------- END OF SCRIPT 1 --------------------

remove(list = ls())
workflow_name <- "otu" # e.g., med, dada2, etc
work_here <- paste0("working_files/ssu/", workflow_name)
share_here <- paste0("share/ssu/", workflow_name)
source(here("assets/scripts", "summarize_objs.R"))
source(here("assets/scripts", "mia_metrics.R"))

#remove(list = ls())
load(here("page_build", "otu_part2.rdata"))
objects()

make.shared(list=final.opti_mcc.list, count=final.count_table, label=0.03)

classify.otu(list=final.opti_mcc.list, count=final.count_table, taxonomy=final.taxonomy, label=0.03)

count.groups(shared=final.opti_mcc.shared)

# pandoc.table(taxonomy_table_16S[10:12, 1:3], emphasize.rownames = FALSE)

# tmp_tax <- read_delim(here(work_here, "final.opti_mcc.0.03.cons.taxonomy"),
#                       delim = "\t")

#tax1 <- tmp_tax
pandoc.table(tax1[11:13, 1:3], emphasize.rownames = FALSE, split.tables = Inf)

# tmp_tax <- data.frame(sapply(tmp_tax,
#                              gsub,
#                              pattern = "\\(\\d+\\)",
#                              replacement = ""))
# tmp_tax <- data.frame(sapply(tmp_tax,
#                              gsub,
#                              pattern = ";$",
#                              replacement = ""))
# tmp_tax <- separate_wider_delim(tmp_tax,
#                               cols = Taxonomy,
#                               delim = ";", names = c(
#                                 "Kingdom", "Phylum",
#                                 "Class", "Order",
#                                 "Family", "Genus"
#                                 )
#                               )
# tmp_tax <- data.frame(sapply(tmp_tax, gsub,
#                            pattern = "^.*_unclassified$",
#                            replacement = ""))
# tmp_tax$Size <- NULL
# tmp_tax <- tibble::column_to_rownames(tmp_tax, "OTU")

#tax2 <- tmp_tax
pandoc.table(tax2[11:13, 1:3], emphasize.rownames = FALSE, split.tables = Inf)

# rank_prefixes <- c(
#   Kingdom = "k",
#   Phylum  = "p",
#   Class   = "c",
#   Order   = "o",
#   Family  = "f",
#   Genus   = "g"
# )
# 
# tmp_tax <- tmp_tax %>%
#   mutate(across(everything(), ~replace_na(., ""))) %>%
#   mutate(across(names(rank_prefixes),
#                 ~ paste0(rank_prefixes[cur_column()], "__", .))) %>%
# tidy_taxonomy()

#tax3 <- tmp_tax
pandoc.table(tax3[11:13, 1:3], emphasize.rownames = FALSE, split.tables = Inf)

pandoc.table(otu_table_16S[1:3, 1:11], emphasize.rownames = FALSE)

# tmp_st <- readr::read_delim(
#   here(work_here, "final.opti_mcc.shared"),  delim = "\t")

# tmp_st$numOtus <- NULL
# tmp_st$label <- NULL
# tmp_st <- tmp_st %>%
#   tidyr::pivot_longer(cols = c(-1), names_to = "tmp") %>%
#   tidyr::pivot_wider(names_from = c(1))
# 
# tmp_st <- tibble::column_to_rownames(tmp_st, "tmp")

# st <- tmp_st
pandoc.table(st[1:5, 1:3], emphasize.rownames = FALSE)

# only need to run this if reading in processed files
# code adds a tab to the beginning of first line
sed '1s/^/\t/' tmp_final.opti_mcc.fixed.shared > final.opti_mcc.fixed.shared

pandoc.table(sample_info_16S[1:3, 1:4], emphasize.rownames = FALSE, split.tables = Inf)

# samdf <- readRDS(here("working_files/ssu/sampledata", "sample_data.rds"))
# 
# samdf <- samdf %>% tibble::column_to_rownames("SampleID")
# samdf$SampleID <- rownames(samdf)
# samdf <- samdf %>% dplyr::relocate(SampleID)
# 
# samdf <- samdf %>%
#   dplyr::filter(
#     stringr::str_starts(SampleID, "Control", negate = TRUE)
#     )

pandoc.table(samdf[1:3, 1:4], emphasize.rownames = FALSE, split.tables = Inf)

# sample_info <- samdf
# tax_tab <- tmp_tax
# otu_tab <- tmp_st

# tmp_me <- microtable$new(sample_table = sample_info,
#                          otu_table = otu_tab,
#                          tax_table = tax_tab)
# tmp_me

get.oturep(column=final.dist, list=final.opti_mcc.list, count=final.count_table, fasta=final.fasta)

seqkit replace -p "\|.*" -r '' final.opti_mcc.0.03.rep.fasta > tmp2.fa
seqkit replace -p ".*\\t" -r '' tmp2.fa > tmp3.fa
seqkit replace -p "-" -r '$1' -s -w 0 tmp3.fa > otu_reps.fasta
rm tmp*

# rep_fasta <- Biostrings::readDNAStringSet(here(work_here, "otu_reps.fasta"))
# tmp_me$rep_fasta <- rep_fasta
# tmp_me$tidy_dataset()
# tmp_me

# me_raw <- microeco::clone(tmp_me)

# threshold <- 1000
# tmp_no_low <- microeco::clone(me_raw)
# tmp_no_low$otu_table <- me_raw$otu_table %>%
#           dplyr::select(where(~ is.numeric(.) && sum(.) >= threshold))
# tmp_no_low$tidy_dataset()
# tmp_no_low

# me_no_low <- microeco::clone(tmp_no_low)

me_no_low

# me_final <- microeco::clone(tmp_no_low)

# ps_final <- file2meco::meco2phyloseq(me_final)

# identical(rownames(me_raw$sample_table), colnames(me_raw$otu_table))
# identical(rownames(me_final$sample_table), colnames(me_final$otu_table))

# objs <- c("me_raw", "me_final")
# pipe_summary <- summarize_objs(objs)
# pipe_summary$me_dataset <- c(
#   "original", "final")
# print(pipe_summary)

knitr::kable(pipe_summary)

# tmp_raw_ps <- file2meco::meco2phyloseq(me_raw)
# tmp_final_ps <- file2meco::meco2phyloseq(me_final)
# tmp_mia_raw <- mia_metrics(tmp_raw_ps)
# tmp_mia_final <- mia_metrics(tmp_final_ps)
# 
# tmp_metrics_final <- rbind(tmp_mia_raw, tmp_mia_final)
# rownames(tmp_metrics_final) <- c("Start", "End")
# 
# # suppose ds_metrics_final has rownames "Start" and "End"
# ds_metrics_tbl <- tmp_metrics_final %>%
#   t() %>%                         # transpose: metrics in rows
#   as.data.frame() %>%
#   tibble::rownames_to_column("Metric") %>%
#   dplyr::rename(Start = Start, End = End)
# 
# ds_metrics_tbl <- ds_metrics_tbl %>%
#   mutate(
#     Start = ifelse(Start %% 1 == 0,
#                    formatC(Start, format = "f", digits = 0),
#                    formatC(Start, format = "f", digits = 3)),
#     End   = ifelse(End %% 1 == 0,
#                    formatC(End, format = "f", digits = 0),
#                    formatC(End, format = "f", digits = 3))
#   )
# 
# ds_metrics_tbl$Metric <- c(
#   "Min. no. of reads",
#   "Max. no. of reads",
#   "Total no. of reads",
#   "Avg. no. of reads",
#   "Median no. of reads",
#   "Total ASVs",
#   "No. of singleton ASVs",
#   "% of singleton ASVs",
#   "Sparsity"
# )

knitr::kable(ds_metrics_tbl, format = "markdown")

# # helper functions
# make_tables <- function(obj, prefix) {
#   rc <- obj$sample_sums() %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("SampleID") %>%
#     dplyr::rename(!!paste0(prefix, "_rc") := 2)
# 
#   asv <- obj$otu_table %>%
#     t() %>%
#     as.data.frame() %>%
#     { rowSums(. > 0) } %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("SampleID") %>%
#     dplyr::rename(!!paste0(prefix, "_otu") := 2)
# 
#   list(rc = rc, asv = asv)
# }

# #tmp_start <- microeco::clone(me_raw)
# #tmp_final <- microeco::clone(me_final)
# # --- input rc file ---
# tmp_rc <- readr::read_delim(here(work_here, "mothur_pipeline_read_changes.txt"))
# 
# # --- generate tables ---
# start_tbls <- make_tables(microeco::clone(me_raw), "start")
# final_tbls <- make_tables(microeco::clone(me_final), "final")
# 
# # --- summary ---
# curate_summary <- tmp_rc %>%
#   dplyr::full_join(start_tbls$rc,  by = "SampleID") %>%
#   dplyr::left_join(start_tbls$asv, by = "SampleID") %>%
#   dplyr::left_join(final_tbls$rc,  by = "SampleID") %>%
#   dplyr::left_join(final_tbls$asv, by = "SampleID")
# 
# removed_samples <- curate_summary %>%
#   dplyr::filter(is.na(final_rc)) %>%
#   dplyr::pull(SampleID)

removed_samples[!grepl("^Control_", removed_samples)]

# workflow_name <- "otu"
# dir.create(here(share_here, paste0(workflow_name, "_curated_data")),
#            recursive = FALSE, showWarnings = TRUE)

# tmp_path <- here(share_here, paste0(workflow_name, "_curated_data/"))
# 
# files <- list(
#   otu = paste0(workflow_name, "_otu_table.txt"),
#   tax = paste0(workflow_name, "_tax_table.txt"),
#   sample = paste0(workflow_name, "_sample_table.txt"),
#   rep = paste0(workflow_name, "_rep.fasta"),
#   me = paste0("me_", workflow_name, ".rds"),
#   ps = paste0("ps_", workflow_name, ".rds"),
#   counts = paste0(workflow_name, "_track_read_counts.txt")
# )

# #---------------OTU Table--------------#
# tmp_otu <- me_final$otu_table %>%
#   tibble::rownames_to_column(paste0(toupper(workflow_name), "_ID"))
# write_delim(tmp_otu, paste0(tmp_path, files$otu), delim = "\t")
# #---------------TAX Table--------------#
# tmp_tax <- me_final$tax_table %>%
#   tibble::rownames_to_column(paste0(toupper(workflow_name), "_ID"))
# write_delim(tmp_tax, paste0(tmp_path, files$tax), delim = "\t")
# #---------------SAMP Table--------------#
# write_delim(me_final$sample_table, paste0(tmp_path, files$sample), delim = "\t")
# #---------------REP FASTA--------------#
# write.fasta(
#   sequences = as.list(as.character(me_final$rep_fasta)),
#   names = names(me_final$rep_fasta),
#   file.out = paste0(tmp_path, files$rep)
# )
# #---------------        --------------#
# saveRDS(me_final, paste0(tmp_path, files$me))
# ps_final <- file2meco::meco2phyloseq(me_final)
# saveRDS(ps_final, paste0(tmp_path, files$ps))
# #---------------read count changes--------------#
# write_delim(curate_summary, paste0(tmp_path, files$counts), delim = "\t")
# file.copy(here("working_files/ssu/sampledata",  "all_metadata.txt"),
#           tmp_path,
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)

# zip::zip(zipfile = here(share_here, paste0(workflow_name, "_curated_data.zip")),
#          files = here(share_here, paste0(workflow_name,"_curated_data")),
#          mode = "cherry-pick")

# file.copy(here(work_here,  "mothur_pipeline_read_changes.txt"),
#           here(share_here),
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)

# dir.create(here(share_here, paste0(workflow_name, "_processing")),
#            recursive = FALSE, showWarnings = TRUE)
# copy_here <- here(share_here, paste0(workflow_name, "_processing/"))

# tmp_sampdata_path <- "working_files/ssu/sampledata"
# fs::dir_copy(path = here(tmp_sampdata_path,  "fastq_rename_results/"),
#              new_path = here(copy_here,  "fastq_rename_results/"),
#              overwrite = TRUE)
# 
# fs::dir_copy(path = here(tmp_sampdata_path,  "fastq_rename_lookup/"),
#              new_path = here(copy_here,  "fastq_rename_lookup/"),
#              overwrite = TRUE)
# file.copy(here(tmp_sampdata_path,  "rename.sh"),
#           here(copy_here),
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)

# fs::dir_copy(path = here(work_here,  paste0(workflow_name, "_hydra_scripts/")),
#              new_path = here(copy_here,  paste0(workflow_name, "_hydra_scripts/")),
#              overwrite = TRUE)
# 
# file.copy(here(work_here,  paste0(workflow_name, "_batchfile_processing/")),
#           copy_here,
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)
# 
# file.copy(here(share_here,  "mothur_pipeline_read_changes.txt"),
#           here(copy_here),
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)
# 
# file.copy(here(tmp_sampdata_path,  "all_metadata.txt"),
#           here(copy_here),
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)

# zip::zip(zipfile = here(share_here, paste0(workflow_name, "_processing.zip")),
#          files = here(share_here, paste0(workflow_name,"_processing")),
#          mode = "cherry-pick")
# dir_delete(here(share_here, paste0(workflow_name,"_processing")))

workflow_name <- "otu"
options(knitr.duplicate.label = "allow")
# Define sources and numbered outputs
sources <- c(
  "include/_OTU_Part1.qmd",
  "index.qmd"
)
# Generate corresponding temporary output files
tmp_outputs <- here("share/ssu", paste0(workflow_name, "_workflow", seq_along(sources), ".R"))
# Final combined output
final_output <- here(share_here, paste0(workflow_name, "_workflow.R"))
# Step 1: Extract code from qmd to temporary R scripts
purrr::walk2(sources, tmp_outputs, ~
  knitr::purl(.x, output = .y, documentation = 0)
)
# Step 2: Concatenate all temporary R scripts into one, with separators
sink(final_output)
purrr::walk2(tmp_outputs, seq_along(tmp_outputs), ~{
  cat(readLines(.x), sep = "\n")
  cat("\n\n",
      paste0("# -------------------- END OF SCRIPT ", .y, " --------------------\n\n"),
      sep = "")
})
sink()
# Step 3: Clean up temporary files
file.remove(tmp_outputs)

# objects()
# gdata::keep(me_final, me_no_low, me_raw, curate_summary,
#             ds_metrics_tbl, pipe_summary, removed_samples,
#             samdf, st, tax1, tax2, tax3,
#             sure = TRUE)
# save.image(here("page_build", "otu_part2.rdata"))


# -------------------- END OF SCRIPT 2 --------------------

