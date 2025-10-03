# project_root <- quarto::find_project_root()
# # Construct the full path to your functions script
# sampdata_path <- file.path(project_root, "working_files/ssu/dada2")

# path <- "BCS_26/"
# head(list.files(path))

c("Control_10_R1_001.trimmed.fastq", "Control_10_R2_001.trimmed.fastq", 
  "Control_11_R1_001.trimmed.fastq", "Control_11_R2_001.trimmed.fastq", 
  "Control_12_R1_001.trimmed.fastq", "Control_12_R2_001.trimmed.fastq")

# fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
# fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
# sample.names <- sapply(strsplit(fnFs, "_R1_"), `[`, 1)
# fnFs <- file.path(path, fnFs)
# fnRs <- file.path(path, fnRs)

# qprofile_fwd <- plotQualityProfile(fnFs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile_rev <- plotQualityProfile(fnRs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile <- grid.arrange(qprofile_fwd, qprofile_rev, nrow = 1)

# ggsave("figures/BCS_26_filt_plot_qscores.png", qprofile, width = 7, height = 3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_26_filt_plot_qscores.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_26_filt_plot_qscores.png")

# filtFs <- file.path(path, "filtered",
#                     paste0(sample.names, "_F_filt.fastq")
#                     )
# filtRs <- file.path(path, "filtered",
#                     paste0(sample.names, "_R_filt.fastq")
#                     )
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#                      truncLen=c(220,180),
#                      maxN = 0, maxEE = 2, truncQ = 2,
#                      rm.phix = TRUE, compress = TRUE,
#                      multithread = 20)

# samptab <- read.table(here(work_here, "tables/BCS_26_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(samptab)[1] <- "SampleID"
# 
# samptab <- samptab %>%
#   select(SampleID, input, filtered)
# 
# samptab$per_reads_kept <- round(samptab$filtered/samptab$input,
#                                   digits = 3)

# errF <- learnErrors(filtFs, multithread = TRUE)

# plotErrors(errF, nominalQ = TRUE)

# p3 <- plotErrors(errF, nominalQ = TRUE)
# ggsave("figures/BCS_26_plot_errorF_1.png", p3, width = 7, height = 5)
# ggsave("figures/BCS_26_plot_errorF_2.png", p3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_26_plot_errorF_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_26_plot_errorF_2.png")

# errR <- learnErrors(filtRs, multithread = TRUE)

# plotErrors(errR, nominalQ = TRUE)

# p4 <- plotErrors(errR, nominalQ = TRUE)
# ggsave("figures/BCS_26_plot_errorR_1.png", p4, width = 7, height = 5)
# ggsave("figures/BCS_26_plot_errorR_2.png", p4)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_26_plot_errorR_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_26_plot_errorR_2.png")

# sam.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
# derepFs <- derepFastq(filtFs)
# names(derepFs) <- sam.names
# derepRs <- derepFastq(filtRs)
# names(derepRs) <- sam.names

# dadaFs <- dada(derepFs, err = errF, pool = "pseudo",
#                multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, pool = "pseudo",
#                multithread = TRUE)

# dadaFs[[1]]

# dadaRs[[1]]

# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# BCS_26 <- makeSequenceTable(mergers)
# dim(BCS_26)

seqtabX <- c(384, 29430)

# table(nchar(getSequences(BCS_26)))

# file.copy(from = paste0(sampdata_path, "/figures/read_length_before_pseudo_BCS_26.png"), to = "figures/")

knitr::include_graphics("include/figures/read_length_before_pseudo_BCS_26.png")

# read_changes <- read.table(here(work_here, "tables/BCS_26_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(read_changes)[1] <- "SampleID"
# 
# read_changes$per_reads_kept <- round(read_changes$merged/read_changes$input,
#                                   digits = 3)
# write_delim(read_changes, (here(work_here, "dada2_read_changes/BCS_26_read_changes.txt")),
#             delim = "\t")

# saveRDS(BCS_26, "BCS_26.rds")


# -------------------- END OF SCRIPT 1 --------------------

# project_root <- quarto::find_project_root()
# # Construct the full path to your functions script
# sampdata_path <- file.path(project_root, "working_files/ssu/dada2")

# path <- "BCS_28/"
# head(list.files(path))

c("Control_31_R1_001.trimmed.fastq", "Control_31_R2_001.trimmed.fastq", 
  "Control_32_R1_001.trimmed.fastq", "Control_32_R2_001.trimmed.fastq", 
  "Control_33_R1_001.trimmed.fastq", "Control_33_R2_001.trimmed.fastq")

# fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
# fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
# sample.names <- sapply(strsplit(fnFs, "_R1_"), `[`, 1)
# fnFs <- file.path(path, fnFs)
# fnRs <- file.path(path, fnRs)

# qprofile_fwd <- plotQualityProfile(fnFs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile_rev <- plotQualityProfile(fnRs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile <- grid.arrange(qprofile_fwd, qprofile_rev, nrow = 1)

# ggsave("figures/BCS_28_filt_plot_qscores.png", qprofile, width = 7, height = 3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_28_filt_plot_qscores.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_28_filt_plot_qscores.png")

# filtFs <- file.path(path, "filtered",
#                     paste0(sample.names, "_F_filt.fastq")
#                     )
# filtRs <- file.path(path, "filtered",
#                     paste0(sample.names, "_R_filt.fastq")
#                     )
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#                      truncLen=c(220,180),
#                      maxN = 0, maxEE = 2, truncQ = 2,
#                      rm.phix = TRUE, compress = TRUE,
#                      multithread = 20)

# samptab <- read.table(here(work_here, "tables/BCS_28_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(samptab)[1] <- "SampleID"
# 
# samptab <- samptab %>%
#   select(SampleID, input, filtered)
# 
# samptab$per_reads_kept <- round(samptab$filtered/samptab$input,
#                                   digits = 3)
# 
# #write_delim(samptab, (here(share_here, "tables/BCS_28_filter_read_changes.txt")),
# #            delim = "\t")

# errF <- learnErrors(filtFs, multithread = TRUE)

# plotErrors(errF, nominalQ = TRUE)

# p3 <- plotErrors(errF, nominalQ = TRUE)
# ggsave("figures/BCS_28_plot_errorF_1.png", p3, width = 7, height = 5)
# ggsave("figures/BCS_28_plot_errorF_2.png", p3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_28_plot_errorF_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_28_plot_errorF_2.png")

# errR <- learnErrors(filtRs, multithread = TRUE)

# plotErrors(errR, nominalQ = TRUE)

# p4 <- plotErrors(errR, nominalQ = TRUE)
# ggsave("figures/BCS_28_plot_errorR_1.png", p4, width = 7, height = 5)
# ggsave("figures/BCS_28_plot_errorR_2.png", p4)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_28_plot_errorR_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_28_plot_errorR_2.png")

# sam.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
# derepFs <- derepFastq(filtFs)
# names(derepFs) <- sam.names
# derepRs <- derepFastq(filtRs)
# names(derepRs) <- sam.names

# dadaFs <- dada(derepFs, err = errF, pool = "pseudo",
#                multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, pool = "pseudo",
#                multithread = TRUE)

# dadaFs[[1]]

# dadaRs[[1]]

# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# BCS_28 <- makeSequenceTable(mergers)
# dim(BCS_28)

seqtabX <- c(192, 9461)

# table(nchar(getSequences(BCS_28)))

# file.copy(from = paste0(sampdata_path, "/figures/read_length_before_pseudo_BCS_28.png"), to = "figures/")

knitr::include_graphics("include/figures/read_length_before_pseudo_BCS_28.png")

# read_changes <- read.table(here(work_here, "tables/BCS_28_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(read_changes)[1] <- "SampleID"
# 
# read_changes$per_reads_kept <- round(read_changes$merged/read_changes$input,
#                                   digits = 3)
# write_delim(read_changes, (here(work_here, "dada2_read_changes/BCS_28_read_changes.txt")),
#             delim = "\t")

# saveRDS(BCS_28, "BCS_28.rds")


# -------------------- END OF SCRIPT 2 --------------------

# project_root <- quarto::find_project_root()
# # Construct the full path to your functions script
# sampdata_path <- file.path(project_root, "working_files/ssu/dada2")

# path <- "BCS_29/"
# head(list.files(path))

c("Control_37_R1_001.trimmed.fastq", "Control_37_R2_001.trimmed.fastq", 
  "Control_38_R1_001.trimmed.fastq", "Control_38_R2_001.trimmed.fastq", 
  "Control_39_R1_001.trimmed.fastq", "Control_39_R2_001.trimmed.fastq")

# fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
# fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
# sample.names <- sapply(strsplit(fnFs, "_R1_"), `[`, 1)
# fnFs <- file.path(path, fnFs)
# fnRs <- file.path(path, fnRs)

# qprofile_fwd <- plotQualityProfile(fnFs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile_rev <- plotQualityProfile(fnRs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile <- grid.arrange(qprofile_fwd, qprofile_rev, nrow = 1)

# ggsave("figures/BCS_29_filt_plot_qscores.png", qprofile, width = 7, height = 3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_29_filt_plot_qscores.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_29_filt_plot_qscores.png")

# filtFs <- file.path(path, "filtered",
#                     paste0(sample.names, "_F_filt.fastq")
#                     )
# filtRs <- file.path(path, "filtered",
#                     paste0(sample.names, "_R_filt.fastq")
#                     )
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#                      truncLen=c(220,180),
#                      maxN = 0, maxEE = 2, truncQ = 2,
#                      rm.phix = TRUE, compress = TRUE,
#                      multithread = 20)

# samptab <- read.table(here(work_here, "tables/BCS_29_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(samptab)[1] <- "SampleID"
# 
# samptab <- samptab %>%
#   select(SampleID, input, filtered)
# 
# samptab$per_reads_kept <- round(samptab$filtered/samptab$input,
#                                   digits = 3)
# 
# #write_delim(samptab, (here(share_here, "tables/BCS_29_filter_read_changes.txt")),
# #            delim = "\t")

# errF <- learnErrors(filtFs, multithread = TRUE)

# plotErrors(errF, nominalQ = TRUE)

# p3 <- plotErrors(errF, nominalQ = TRUE)
# ggsave("figures/BCS_29_plot_errorF_1.png", p3, width = 7, height = 5)
# ggsave("figures/BCS_29_plot_errorF_2.png", p3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_29_plot_errorF_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_29_plot_errorF_2.png")

# errR <- learnErrors(filtRs, multithread = TRUE)

# plotErrors(errR, nominalQ = TRUE)

# p4 <- plotErrors(errR, nominalQ = TRUE)
# ggsave("figures/BCS_29_plot_errorR_1.png", p4, width = 7, height = 5)
# ggsave("figures/BCS_29_plot_errorR_2.png", p4)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_29_plot_errorR_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_29_plot_errorR_2.png")

# sam.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
# derepFs <- derepFastq(filtFs)
# names(derepFs) <- sam.names
# derepRs <- derepFastq(filtRs)
# names(derepRs) <- sam.names

# dadaFs <- dada(derepFs, err = errF, pool = "pseudo",
#                multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, pool = "pseudo",
#                multithread = TRUE)

# dadaFs[[1]]

# dadaRs[[1]]

# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# BCS_29 <- makeSequenceTable(mergers)
# dim(BCS_29)

seqtabX <- c(384, 21105)

# table(nchar(getSequences(BCS_29)))

# file.copy(from = paste0(sampdata_path, "/figures/read_length_before_pseudo_BCS_29.png"), to = "figures/")

knitr::include_graphics("include/figures/read_length_before_pseudo_BCS_29.png")

# read_changes <- read.table(here(work_here, "tables/BCS_29_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(read_changes)[1] <- "SampleID"
# 
# read_changes$per_reads_kept <- round(read_changes$merged/read_changes$input,
#                                   digits = 3)
# write_delim(read_changes, (here(work_here, "dada2_read_changes/BCS_29_read_changes.txt")),
#             delim = "\t")

# saveRDS(BCS_29, "BCS_29.rds")


# -------------------- END OF SCRIPT 3 --------------------

# project_root <- quarto::find_project_root()
# # Construct the full path to your functions script
# sampdata_path <- file.path(project_root, "working_files/ssu/dada2")

# path <- "BCS_30/"
# head(list.files(path))

c("Control_49_R1_001.trimmed.fastq", "Control_49_R2_001.trimmed.fastq", 
  "Control_50_R1_001.trimmed.fastq", "Control_50_R2_001.trimmed.fastq", 
  "Control_51_R1_001.trimmed.fastq", "Control_51_R2_001.trimmed.fastq")

# fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
# fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
# sample.names <- sapply(strsplit(fnFs, "_R1_"), `[`, 1)
# fnFs <- file.path(path, fnFs)
# fnRs <- file.path(path, fnRs)

# qprofile_fwd <- plotQualityProfile(fnFs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile_rev <- plotQualityProfile(fnRs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile <- grid.arrange(qprofile_fwd, qprofile_rev, nrow = 1)

# ggsave("figures/BCS_30_filt_plot_qscores.png", qprofile, width = 7, height = 3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_30_filt_plot_qscores.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_30_filt_plot_qscores.png")

# filtFs <- file.path(path, "filtered",
#                     paste0(sample.names, "_F_filt.fastq")
#                     )
# filtRs <- file.path(path, "filtered",
#                     paste0(sample.names, "_R_filt.fastq")
#                     )
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#                      truncLen=c(220,180),
#                      maxN = 0, maxEE = 2, truncQ = 2,
#                      rm.phix = TRUE, compress = TRUE,
#                      multithread = 20)

# samptab <- read.table(here(work_here, "tables/BCS_30_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(samptab)[1] <- "SampleID"
# 
# samptab <- samptab %>%
#   select(SampleID, input, filtered)
# 
# samptab$per_reads_kept <- round(samptab$filtered/samptab$input,
#                                   digits = 3)
# 
# #write_delim(samptab, (here(share_here, "tables/BCS_30_filter_read_changes.txt")),
# #            delim = "\t")

# errF <- learnErrors(filtFs, multithread = TRUE)

# plotErrors(errF, nominalQ = TRUE)

# p3 <- plotErrors(errF, nominalQ = TRUE)
# ggsave("figures/BCS_30_plot_errorF_1.png", p3, width = 7, height = 5)
# ggsave("figures/BCS_30_plot_errorF_2.png", p3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_30_plot_errorF_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_30_plot_errorF_2.png")

# errR <- learnErrors(filtRs, multithread = TRUE)

# plotErrors(errR, nominalQ = TRUE)

# p4 <- plotErrors(errR, nominalQ = TRUE)
# ggsave("figures/BCS_30_plot_errorR_1.png", p4, width = 7, height = 5)
# ggsave("figures/BCS_30_plot_errorR_2.png", p4)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_30_plot_errorR_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_30_plot_errorR_2.png")

# sam.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
# derepFs <- derepFastq(filtFs)
# names(derepFs) <- sam.names
# derepRs <- derepFastq(filtRs)
# names(derepRs) <- sam.names

# dadaFs <- dada(derepFs, err = errF, pool = "pseudo",
#                multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, pool = "pseudo",
#                multithread = TRUE)

# dadaFs[[1]]

# dadaRs[[1]]

# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# BCS_30 <- makeSequenceTable(mergers)
# dim(BCS_30)

seqtabX <- c(380, 36401)

# table(nchar(getSequences(BCS_30)))

# file.copy(from = paste0(sampdata_path, "/figures/read_length_before_pseudo_BCS_30.png"), to = "figures/")

knitr::include_graphics("include/figures/read_length_before_pseudo_BCS_30.png")

# read_changes <- read.table(here(work_here, "tables/BCS_30_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(read_changes)[1] <- "SampleID"
# 
# read_changes$per_reads_kept <- round(read_changes$merged/read_changes$input,
#                                   digits = 3)
# write_delim(read_changes, (here(work_here, "dada2_read_changes/BCS_30_read_changes.txt")),
#             delim = "\t")

# saveRDS(BCS_30, "BCS_30.rds")


# -------------------- END OF SCRIPT 4 --------------------

# project_root <- quarto::find_project_root()
# # Construct the full path to your functions script
# sampdata_path <- file.path(project_root, "working_files/ssu/dada2")

# path <- "BCS_34/"
# head(list.files(path))

c("Control_25_R1_001.trimmed.fastq", "Control_25_R2_001.trimmed.fastq", 
  "Control_26_R1_001.trimmed.fastq", "Control_26_R2_001.trimmed.fastq", 
  "Control_27_R1_001.trimmed.fastq", "Control_27_R2_001.trimmed.fastq")

# fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
# fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
# sample.names <- sapply(strsplit(fnFs, "_R1_"), `[`, 1)
# fnFs <- file.path(path, fnFs)
# fnRs <- file.path(path, fnRs)

# qprofile_fwd <- plotQualityProfile(fnFs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile_rev <- plotQualityProfile(fnRs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile <- grid.arrange(qprofile_fwd, qprofile_rev, nrow = 1)

# ggsave("figures/BCS_34_filt_plot_qscores.png", qprofile, width = 7, height = 3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_34_filt_plot_qscores.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_34_filt_plot_qscores.png")

# filtFs <- file.path(path, "filtered",
#                     paste0(sample.names, "_F_filt.fastq")
#                     )
# filtRs <- file.path(path, "filtered",
#                     paste0(sample.names, "_R_filt.fastq")
#                     )
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#                      truncLen=c(220,180),
#                      maxN = 0, maxEE = 2, truncQ = 2,
#                      rm.phix = TRUE, compress = TRUE,
#                      multithread = 20)

# samptab <- read.table(here(work_here, "tables/BCS_34_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(samptab)[1] <- "SampleID"
# 
# samptab <- samptab %>%
#   select(SampleID, input, filtered)
# 
# samptab$per_reads_kept <- round(samptab$filtered/samptab$input,
#                                   digits = 3)
# 
# #write_delim(samptab, (here(share_here, "tables/BCS_34_filter_read_changes.txt")),
# #            delim = "\t")

# errF <- learnErrors(filtFs, multithread = TRUE)

# plotErrors(errF, nominalQ = TRUE)

# p3 <- plotErrors(errF, nominalQ = TRUE)
# ggsave("figures/BCS_34_plot_errorF_1.png", p3, width = 7, height = 5)
# ggsave("figures/BCS_34_plot_errorF_2.png", p3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_34_plot_errorF_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_34_plot_errorF_2.png")

# errR <- learnErrors(filtRs, multithread = TRUE)

# plotErrors(errR, nominalQ = TRUE)

# p4 <- plotErrors(errR, nominalQ = TRUE)
# ggsave("figures/BCS_34_plot_errorR_1.png", p4, width = 7, height = 5)
# ggsave("figures/BCS_34_plot_errorR_2.png", p4)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_34_plot_errorR_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_34_plot_errorR_2.png")

# sam.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
# derepFs <- derepFastq(filtFs)
# names(derepFs) <- sam.names
# derepRs <- derepFastq(filtRs)
# names(derepRs) <- sam.names

# dadaFs <- dada(derepFs, err = errF, pool = "pseudo",
#                multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, pool = "pseudo",
#                multithread = TRUE)

# dadaFs[[1]]

# dadaRs[[1]]

# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# BCS_34 <- makeSequenceTable(mergers)
# dim(BCS_34)

seqtabX <- c(190, 18373)

# table(nchar(getSequences(BCS_34)))

# file.copy(from = paste0(sampdata_path, "/figures/read_length_before_pseudo_BCS_34.png"), to = "figures/")

knitr::include_graphics("include/figures/read_length_before_pseudo_BCS_34.png")

# read_changes <- read.table(here(work_here, "tables/BCS_34_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(read_changes)[1] <- "SampleID"
# 
# read_changes$per_reads_kept <- round(read_changes$merged/read_changes$input,
#                                   digits = 3)
# write_delim(read_changes, (here(work_here, "dada2_read_changes/BCS_34_read_changes.txt")),
#             delim = "\t")

# saveRDS(BCS_34, "BCS_34.rds")


# -------------------- END OF SCRIPT 5 --------------------

# project_root <- quarto::find_project_root()
# # Construct the full path to your functions script
# sampdata_path <- file.path(project_root, "working_files/ssu/dada2")

# path <- "BCS_35/"
# head(list.files(path))

c("Control_13_R1_001.trimmed.fastq","Control_13_R2_001.trimmed.fastq", 
  "Control_14_R1_001.trimmed.fastq", "Control_14_R2_001.trimmed.fastq", 
  "Control_15_R1_001.trimmed.fastq", "Control_15_R2_001.trimmed.fastq")

# fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
# fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
# sample.names <- sapply(strsplit(fnFs, "_R1_"), `[`, 1)
# fnFs <- file.path(path, fnFs)
# fnRs <- file.path(path, fnRs)

# qprofile_fwd <- plotQualityProfile(fnFs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile_rev <- plotQualityProfile(fnRs[1:x],
#                                    aggregate = TRUE,
#                                    n = 20000)
# qprofile <- grid.arrange(qprofile_fwd, qprofile_rev, nrow = 1)

# ggsave("figures/BCS_35_filt_plot_qscores.png", qprofile, width = 7, height = 3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_35_filt_plot_qscores.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_35_filt_plot_qscores.png")

# filtFs <- file.path(path, "filtered",
#                     paste0(sample.names, "_F_filt.fastq")
#                     )
# filtRs <- file.path(path, "filtered",
#                     paste0(sample.names, "_R_filt.fastq")
#                     )
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
#                      truncLen=c(220,180),
#                      maxN = 0, maxEE = 2, truncQ = 2,
#                      rm.phix = TRUE, compress = TRUE,
#                      multithread = 20)

# samptab <- read.table(here(work_here, "tables/BCS_35_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(samptab)[1] <- "SampleID"
# 
# samptab <- samptab %>%
#   select(SampleID, input, filtered)
# 
# samptab$per_reads_kept <- round(samptab$filtered/samptab$input,
#                                   digits = 3)
# #write_delim(samptab, (here(share_here, "tables/BCS_35_filter_read_changes.txt")),
# #            delim = "\t")

# errF <- learnErrors(filtFs, multithread = TRUE)

# plotErrors(errF, nominalQ = TRUE)

# p3 <- plotErrors(errF, nominalQ = TRUE)
# ggsave("figures/BCS_35_plot_errorF_1.png", p3, width = 7, height = 5)
# ggsave("figures/BCS_35_plot_errorF_2.png", p3)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_35_plot_errorF_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_35_plot_errorF_2.png")

# errR <- learnErrors(filtRs, multithread = TRUE)

# plotErrors(errR, nominalQ = TRUE)

# p4 <- plotErrors(errR, nominalQ = TRUE)
# ggsave("figures/BCS_35_plot_errorR_1.png", p4, width = 7, height = 5)
# ggsave("figures/BCS_35_plot_errorR_2.png", p4)

# file.copy(from = paste0(sampdata_path, "/figures/BCS_35_plot_errorR_2.png"), to = "figures/")

knitr::include_graphics("include/figures/BCS_35_plot_errorR_2.png")

# sam.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
# derepFs <- derepFastq(filtFs)
# names(derepFs) <- sam.names
# derepRs <- derepFastq(filtRs)
# names(derepRs) <- sam.names

# dadaFs <- dada(derepFs, err = errF, pool = "pseudo",
#                multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, pool = "pseudo",
#                multithread = TRUE)

# dadaFs[[1]]

# dadaRs[[1]]

# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# BCS_35 <- makeSequenceTable(mergers)
# dim(BCS_35)

seqtabX <- c(379, 23367)

# table(nchar(getSequences(BCS_35)))

# file.copy(from = paste0(sampdata_path, "/figures/read_length_before_pseudo_BCS_35.png"), to = "figures/")

knitr::include_graphics("include/figures/read_length_before_pseudo_BCS_35.png")

# read_changes <- read.table(here(work_here, "tables/BCS_35_read_changes.txt"),
#                        header = TRUE, sep = "\t")
# names(read_changes)[1] <- "SampleID"
# 
# read_changes$per_reads_kept <- round(read_changes$merged/read_changes$input,
#                                   digits = 3)
# write_delim(read_changes, (here(work_here, "dada2_read_changes/BCS_35_read_changes.txt")),
#             delim = "\t")

# saveRDS(BCS_35, "BCS_35.rds")


# -------------------- END OF SCRIPT 6 --------------------

# project_root <- quarto::find_project_root()
# # Construct the full path to your functions script
# sampdata_path <- file.path(project_root, "working_files/ssu/dada2")

# BCS_26 <- readRDS("`BCS_26.rds")
# BCS_28 <- readRDS("`BCS_28.rds")
# BCS_29 <- readRDS("`BCS_29.rds")
# BCS_30 <- readRDS("`BCS_30.rds")
# BCS_34 <- readRDS("`BCS_34.rds")
# BCS_35 <- readRDS("`BCS_35.rds")

# seqtab.merge <- mergeSequenceTables(BCS_26, BCS_28, BCS_29,
#                                     BCS_30, BCS_34, BCS_35)
# dim(seqtab.merge)

seqtabZ <- c(1909, 96680)
seqtabZ

# table(nchar(getSequences(seqtab.merge)))

# read_length_all <-  data.frame(nchar(getSequences(seqtab.merge)))
# colnames(read_length_all) <- "length"
# plot_all <- qplot(length, data = read_length_all, geom = "histogram",
#                   binwidth = 1, xlab = "read length",
#                   ylab = "total variants", xlim = c(200,400))

# ggsave("read_length_before_collapse.png", plot_all, width = 7, height = 3)

# file.copy(from = paste0(sampdata_path, "/figures/read_length_before_collapse.png"), to = "figures/")

knitr::include_graphics("include/figures/read_length_before_collapse.png")

# seqtab.trim <- seqtab.merge[,nchar(colnames(seqtab.merge)) %in%
#                               seq(252, 254)]
# dim(seqtab.trim)

seqtabY <- c(1909, 91249)
seqtabY

# table(nchar(getSequences(seqtab.trim)))

# seqtab.trim.nochim.consensus <-
#   removeBimeraDenovo(seqtab.trim,
#                      method = "consensus",
#                      multithread = 20,
#                      verbose = TRUE)
# dim(seqtab.trim.nochim.consensus)

seqtab3 <- c(1909, 72851)
seqtab3

# sum(seqtab.nochim)/sum(seqtab.2)

chim_rem <- 0.9669996

# getN <- function(x) sum(getUniques(x))
# track <- cbind(rowSums(seqtab),
#                rowSums(seqtab.trim),
#                rowSums(seqtab.trim.nochim.pool),
#                rowSums(seqtab.trim.nochim.consensus))
# 
# colnames(track) <- c("merged", "trim",
#                      "chimera_pool",
#                      "chimera_concensus")

# tmp_listfile <- list.files(path = here(work_here, "tables"),
#                            pattern = "_read_changes.txt",
#                            full.names = TRUE, recursive = FALSE)
# 
# read_file <- function(filename) {
#   dat <- read.table(filename, header = TRUE, sep = "\t")
#   names(dat)[1] <- "SampleID"
#   return(dat)
# }
# tmp_dat <- do.call(rbind, lapply(tmp_listfile, read_file))
# tmp_tab <- read.table(here(work_here, "tables/3.chimera_read_changes_pipeline.txt"),
#                        header = TRUE, sep = "\t")
# names(tmp_tab)[1] <- "SampleID"
# 
# tmp_full_tab <- dplyr::left_join(tmp_dat, tmp_tab,
#                                  by = c("SampleID" = "SampleID",
#                                         "merged" = "merged"))
# tmp_full_tab[8] <- NULL
# names(tmp_full_tab)[8] <- "nochim"
# write_delim(tmp_full_tab,
#             here(work_here, "dada2_read_changes/all_sample_dada2_read_changes.txt"),
#             delim = "\t")

wget https://manichanh.vhir.org/gsrdb/GSR-DB_V4_cluster-1.tar.gz
tar -xvzf GSR-DB_V4_cluster-1.tar.gz

head GSR-DB_V4_cluster-1_taxa.txt

head GSR-DB_V4_cluster-1_seqs.fasta

conda activate seqkit
seqkit replace -w 0  -p "(.+)" -r '{kv}' -k GSR-DB_V4_cluster-1_taxa.txt GSR-DB_V4_cluster-1_seqs.fasta > tmp_1.fa

seqkit replace -w 0  -p " s__.*" -r ''  tmp_1.fa > tmp_2.fa
seqkit replace -w 0  -p "\s" -r ''  tmp_2.fa > tmp_3.fa
seqkit replace -w 0  -p "\w__" -r ''  tmp_3.fa > gsrdb_dada2.fa
rm tmp_*

# seqtab.consensus <- seqtab.trim.nochim.consensus
# tax_gsrdb.consensus <-
#   assignTaxonomy(seqtab.consensus,
#                  "TAXONOMY_FILES/gsrdb_dada2.fa",
#                  multithread = TRUE,
#                  verbose = TRUE)
# saveRDS(tax_gsrdb.consensus, "4.tax_gsrdb.consensus.rds")


# -------------------- END OF SCRIPT 7 --------------------

remove(list = ls())
load(here("page_build", "dada2_part2.rdata"))
workflow_name <- "dada2" # e.g., med, dada2, etc
work_here <- paste0("working_files/ssu/", workflow_name)
share_here <- paste0("share/ssu/", workflow_name)
source(here("assets/scripts", "summarize_objs.R"))
source(here("assets/scripts", "mia_metrics.R"))

# #!/usr/bin/env Rscript
# set.seed(919191)
# pacman::p_load(tidyverse, gridExtra, grid, phyloseq,
#                formatR, gdata, ff, decontam, dada2,
#                ShortRead, Biostrings, DECIPHER,
#                install = FALSE, update = FALSE)

# seqtab <- readRDS(here(work_here, "rdata/3.seqtab.trim.nochim.consensus.rds"))
# tax_gsrdb <- readRDS(here(work_here, "rdata/4.tax_gsrdb.consensus.rds"))

# tmp_tab1 <- readRDS(here("working_files/ssu", "sampledata/sample_data.rds"))
# tmp_tab2 <- read.table(here(work_here,
#     "dada2_pipeline_read_changes.txt"),
#     header = TRUE, sep = "\t"
# )
# tmp_tab2[3:7] <- NULL
# tmp_tab1 <- arrange(tmp_tab1, SampleID, .by_group = FALSE)
# tmp_tab2 <- arrange(tmp_tab2, SampleID, .by_group = FALSE)
# identical(tmp_tab1$SampleID, tmp_tab2$SampleID)

# tmp_tab3 <- data.frame(row.names(seqtab))
# colnames(tmp_tab3) <- "SampleID"
# 
# tmp_seqtab <- data.frame(seqtab)
# tmp_seqtab <- tibble::rownames_to_column(tmp_seqtab, var = "SampleID")
# 
# identical(tmp_tab3$SampleID, tmp_seqtab$SampleID)

# tmp_seqtab <- tmp_seqtab %>%
#     mutate(count = rowSums(. != 0))
# tmp_tab3$no_asvs <- tmp_seqtab$count

# tmp_tab4 <- dplyr::right_join(tmp_tab1, tmp_tab2, by = "SampleID")
# 
# tmp_tab2 <- tmp_tab2[order(tmp_tab2$SampleID),]
# tmp_tab4 <- tmp_tab4[order(tmp_tab4$SampleID),]
# 
# identical(tmp_tab2$input, tmp_tab4$input)
# identical(tmp_tab2$SampleID, tmp_tab4$SampleID)
# 
# tmp_tab3 <- tmp_tab3[order(tmp_tab3$SampleID),]
# 
# tmp_tab5 <- dplyr::left_join(tmp_tab4, tmp_tab3, by = "SampleID")
# tmp_tab5$per_reads_kept <- round(tmp_tab5$nochim/tmp_tab5$input, digits = 3)
# samptab <- tmp_tab5

# samptab <- samptab %>%
#     dplyr::relocate(
#         c(input, nochim, per_reads_kept, no_asvs),
#         .after = "SampleID"
#     )
# samptab <- samptab %>%
#     dplyr::relocate(per_reads_kept, .after = "nochim")
# 
# rm(list = ls(pattern = "tmp_"))

pandoc.table(taxonomy_table_16S[10:12, 1:3], emphasize.rownames = FALSE)

pandoc.table(data.frame(tax_gsrdb)[1:3, 1:3], emphasize.rownames = FALSE)

# tmp_tax <- data.frame(tax_gsrdb)
# # adding unique ASV names
# row.names(tmp_tax) <- paste0("ASV", seq(nrow(tmp_tax)))

# tax.head1 <- tmp_tax[1:3, 1:3]

pandoc.table(tax.head1, emphasize.rownames = FALSE)

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

# tax.head2 <- tmp_tax[1:3, 1:3]

pandoc.table(tax.head2, emphasize.rownames = FALSE)

pandoc.table(otu_table_16S[1:3, 1:11], emphasize.rownames = FALSE)

# tmp_st <- data.frame(seqtab)
# identical(colnames(tmp_st), row.names(tax_gsrdb))
# names(tmp_st) <- row.names(tmp_tax)
# 
# tmp_st <-  tmp_st %>% tibble::rownames_to_column()
# 
# tmp_st <- tmp_st %>%
#   tidyr::pivot_longer(cols = c(-1), names_to = "tmp") %>%
#   tidyr::pivot_wider(names_from = c(1))
# tmp_st <- tibble::column_to_rownames(tmp_st, "tmp")

# st.head <- tmp_st[1:3, 61:63]

pandoc.table(st.head, emphasize.rownames = FALSE)

pandoc.table(sample_info_16S[1:3,], emphasize.rownames = FALSE)

# samdf <- readRDS(here("share/ssu/sampledata/", "sample_data.rds"))
# samdf <- samdf %>% tibble::column_to_rownames("SampleID")
# samdf$SampleID <- rownames(samdf)
# samdf <- samdf %>% relocate(SampleID)

pandoc.table(samdf[61:63, 1:5], emphasize.rownames = FALSE)

# sample_info <- samdf
# tax_tab <- tmp_tax
# otu_tab <- tmp_st

# tmp_me <- microtable$new(sample_table = sample_info,
#                          otu_table = otu_tab,
#                          tax_table = tax_tab)
# tmp_me

# tmp_seq <- data.frame(row.names(data.frame(tax_gsrdb)) )
# tmp_names <- data.frame(row.names(tax_tab))
# tmp_fa <- cbind(tmp_names, tmp_seq)
# colnames(tmp_fa) <- c("ASV_ID", "ASV_SEQ")
# tmp_fa$ASV_ID <- sub("^", ">", tmp_fa$ASV_ID)
# 
# write.table(tmp_fa, here(work_here, "rep_seq.fasta"),
#             sep = "\n", col.names = FALSE, row.names = FALSE,
#             quote = FALSE, fileEncoding = "UTF-8")
# rep_fasta <- Biostrings::readDNAStringSet(here(work_here, "rep_seq.fasta"))
# tmp_me$rep_fasta <- rep_fasta
# tmp_me$tidy_dataset()

# me_raw <- microeco::clone(tmp_me)

# tmp_no_na <- microeco::clone(tmp_me)
# tmp_no_na$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
# tmp_no_na$tidy_dataset()

# me_no_na <- microeco::clone(tmp_no_na)

me_no_na

# tmp_no_cont <- microeco::clone(tmp_no_na)
# tmp_no_cont$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# tmp_no_cont$tidy_dataset()
# tmp_no_cont

# me_no_cont <- microeco::clone(tmp_no_cont)

me_no_cont

# tmp_nc <- microeco::clone(tmp_no_cont)
# tmp_nc$sample_table <- subset(tmp_nc$sample_table, TAXON == "Control")
# tmp_nc$tidy_dataset()
# tmp_nc

# me_nc <- microeco::clone(tmp_nc)

me_nc

# nc_asvs <- row.names(tmp_nc$tax_table)

head(nc_asvs, n = 20)

# tmp_rem_nc <- microeco::clone(tmp_no_cont)
# tmp_rem_nc_df <- tmp_rem_nc$otu_table
# tmp_rem_nc_df <- tmp_rem_nc_df %>%
#                  dplyr::filter(row.names(tmp_rem_nc_df) %in% nc_asvs)
# tmp_rem_nc_df <- tmp_rem_nc_df %>% tibble::rownames_to_column("ASV_ID")

# #-------provide a string unique to NC samples--------------#
# nc_name <- "Control_"
# #----------------------------------------------------------#
# tmp_rem_nc_df <- tmp_rem_nc_df  %>%
#   dplyr::mutate(total_reads_NC = rowSums(dplyr::select(., contains(nc_name))),
#          .after = "ASV_ID")
# tmp_rem_nc_df <- dplyr::select(tmp_rem_nc_df, -contains(nc_name))
# tmp_rem_nc_df <- tmp_rem_nc_df %>%
#   dplyr::mutate(total_reads_samps = rowSums(.[3:ncol(tmp_rem_nc_df)]),
#                 .after = "total_reads_NC")
# tmp_rem_nc_df[, 4:ncol(tmp_rem_nc_df)] <- list(NULL)
# tmp_rem_nc_df <- tmp_rem_nc_df %>%
#   dplyr::mutate(perc_in_neg = 100*(
#     total_reads_NC / (total_reads_NC + total_reads_samps)),
#                 .after = "total_reads_samps")

# tmp_rem_nc_df$perc_in_neg <- round(tmp_rem_nc_df$perc_in_neg, digits = 6)
# 
# tmp_1 <- data.frame(rowSums(tmp_rem_nc$otu_table != 0)) %>%
#                    tibble::rownames_to_column("ASV_ID") %>%
#                    dplyr::rename("total_samples" = 2)
# 
# tmp_2 <- dplyr::select(tmp_rem_nc$otu_table, contains(nc_name))
# tmp_2$num_samp_nc <- rowSums(tmp_2 != 0)
# tmp_2 <- dplyr::select(tmp_2, contains("num_samp_nc")) %>%
#                       tibble::rownames_to_column("ASV_ID")
# 
# tmp_3 <- dplyr::select(tmp_rem_nc$otu_table, -contains(nc_name))
# tmp_3$num_samp_no_nc <- rowSums(tmp_3 != 0)
# tmp_3 <- dplyr::select(tmp_3, contains("num_samp_no_nc")) %>%
#                       tibble::rownames_to_column("ASV_ID")
# 
# tmp_rem_nc_df <- dplyr::left_join(tmp_rem_nc_df, tmp_1) %>%
#                  dplyr::left_join(., tmp_2) %>%
#                  dplyr::left_join(., tmp_3)
# 
# tmp_rem_nc_df <- tmp_rem_nc_df %>%
#   dplyr::mutate(perc_in_neg_samp = 100*( num_samp_nc / (num_samp_nc + num_samp_no_nc)),
#                 .after = "num_samp_no_nc")

# nc_check <- tmp_rem_nc_df

# nc_remove <- nc_check %>%
#   dplyr::filter(perc_in_neg > 10 | perc_in_neg_samp > 10)

# nc_remain <- dplyr::anti_join(nc_check, nc_remove)
# 
# rem_nc_reads <- sum(nc_remove$total_reads_NC)
# rem_sam_reads <- sum(nc_remove$total_reads_samps)
# per_reads_rem <- round(100*( rem_nc_reads / (rem_nc_reads + rem_sam_reads)),
#                        digits = 3)
# 
# ret_nc_reads <- sum(nc_remain$total_reads_NC)
# ret_sam_reads <- sum(nc_remain$total_reads_samps)
# per_reads_ret <- round(100*( ret_nc_reads / (ret_nc_reads + ret_sam_reads)),
#                        digits = 3)

# tmp_no_nc <- microeco::clone(tmp_no_cont)
# 
# tmp_rem_asv <- as.factor(nc_remove$ASV_ID)
# tmp_no_nc$otu_table <- tmp_rem_nc$otu_table %>%
#   filter(!row.names(tmp_no_nc$otu_table) %in% tmp_rem_asv)
# tmp_no_nc$tidy_dataset()
# 
# tmp_no_nc$sample_table <- subset(tmp_no_nc$sample_table,
#                                  TAXON != "Control_")
# tmp_no_nc$tidy_dataset()
# tmp_no_nc

# me_no_nc <- microeco::clone(tmp_no_nc)

me_no_nc

# tmp_no_low <- microeco::clone(tmp_no_nc)
# tmp_no_low$otu_table <- tmp_no_nc$otu_table %>%
#           dplyr::select(where(~ is.numeric(.) && sum(.) >= 1000))
# tmp_no_low$tidy_dataset()
# tmp_no_low

# me_no_low <- microeco::clone(tmp_no_low)

me_no_low

# me_final <- microeco::clone(tmp_no_low)

# ps_final <- file2meco::meco2phyloseq(me_final)

# identical(rownames(me_raw$sample_table), colnames(me_raw$otu_table))
# identical(rownames(me_final$sample_table), colnames(me_final$otu_table))

# objs <- c("me_raw", "me_no_na",
#           "me_no_cont", "me_no_nc",
#           "me_no_low", "me_final")
# pipe_summary <- summarize_objs(objs)
# pipe_summary$me_dataset <- c(
#   "original", "no NA kingdoms",
#   "no contaminants", "no negative controls",
#   "no low count samps", "final")
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
#     dplyr::rename(!!paste0(prefix, "_asv") := 2)
# 
#   list(rc = rc, asv = asv)
# }

# # --- input rc file ---
# tmp_rc <- readr::read_delim(here(work_here, "dada2_pipeline_read_changes.txt"))
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

# workflow_name <- "dada2"
# dir.create(here(share_here, paste0(workflow_name, "_curated_data")),
#            recursive = FALSE, showWarnings = TRUE)

# write_delim(nc_check, here(share_here, "asv_in_nc_samples.txt"),
#     delim = "\t")

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
#   tibble::rownames_to_column("ASV_ID")
# write_delim(tmp_otu, paste0(tmp_path, files$otu), delim = "\t")
# #---------------TAX Table--------------#
# tmp_tax <- me_final$tax_table %>%
#   tibble::rownames_to_column("ASV_ID")
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

# file.copy(here(work_here,  "dada2_pipeline_read_changes.txt"),
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
# 
# file.copy(here(tmp_sampdata_path,  "rename.sh"),
#           here(copy_here),
#           overwrite = TRUE, recursive = FALSE,
#           copy.mode = TRUE, copy.date = FALSE)

# fs::dir_copy(path = here(work_here,  paste0(workflow_name, "_hydra_scripts/")),
#              new_path = here(copy_here,  paste0(workflow_name, "_hydra_scripts/")),
#              overwrite = TRUE)
# 
# fs::dir_copy(path = here(work_here,  paste0(workflow_name, "_scripts/")),
#              new_path = here(copy_here,  paste0(workflow_name, "_scripts/")),
#              overwrite = TRUE)
# 
# fs::dir_copy(path = here(work_here,  paste0(workflow_name, "_read_changes/")),
#              new_path = here(copy_here,  paste0(workflow_name, "_read_changes/")),
#              overwrite = TRUE)
# 
# fs::dir_copy(path = here(work_here,  "dada_results/"),
#              new_path = here(copy_here,  "dada_results/"),
#              overwrite = TRUE)
# 
# file.copy(here(share_here,  "dada2_pipeline_read_changes.txt"),
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

workflow_name <- "dada2"
options(knitr.duplicate.label = "allow")
# Define sources and numbered outputs
sources <- c(
  "include/_BCS_26.qmd",
  "include/_BCS_28.qmd",
  "include/_BCS_29.qmd",
  "include/_BCS_30.qmd",
  "include/_BCS_34.qmd",
  "include/_BCS_35.qmd",
  "include/_MERGE_RUNS.qmd",
  "index.qmd"
)
# Generate corresponding temporary output files
tmp_outputs <- here(share_here, paste0(workflow_name, "_workflow", seq_along(sources), ".R"))
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
# gdata::keep(curate_summary, ds_metrics_tbl, me_final, me_nc,
#             me_no_cont, me_no_low, me_no_na, me_no_nc,
#             me_raw, nc_asvs, nc_check, nc_name, nc_remain, nc_remove,
#             per_reads_rem, per_reads_ret, pipe_summary, rem_nc_reads,
#             rem_sam_reads, removed_samples, ret_nc_reads, ret_sam_reads,
#             samdf, sample_info, samptab, seqtab, share_here, st.head,
#             tax.head1, tax.head2, seqtab, tax_gsrdb,
#             sure = TRUE)
# save.image(here("page_build", "dada2_part2.rdata"))


# -------------------- END OF SCRIPT 8 --------------------

