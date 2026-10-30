############################################################
# Physalia 2026: Causal AI Methods for Computational Biology
# Build script for the course dataset: iHMP_causal
#
# Inputs (folder raw/):
#   HMP2_metadata.tsv, HMP2_taxonomy.tsv
#       HMP2 stool metagenome metadata and species profiles as
#       distributed with MaAsLin 3 (github.com/biobakery/maaslin3)
#   hmp2_metadata_2018-08-20.csv   full IBDMDB sample metadata
#   HMP2_metabolites.tsv           stool metabolomics (HMP2 MBX)
#   host_tx_counts.pcl.tsv         biopsy host RNA-seq counts (HTX)
#   FUT2.vcf                       host exome genotypes at FUT2
#   HMP.and1kG.MDS.mds             ancestry MDS (exome + 1000 Genomes)
#
# Study: Lloyd-Price J et al. (2019) Multi-omics of the gut
#   microbial ecosystem in inflammatory bowel diseases. Nature
#   569: 655-662. Data portal: https://ibdmdb.org
#
# Output (folder out/):
#   iHMP_causal.rds  named list of the tables below
#   visits.csv, species.csv, metabolites.csv, biopsy.csv,
#   genotype.csv, dictionary.csv
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(vegan)
})

dir.create("out", showWarnings = FALSE)
big <- read.csv("raw/hmp2_metadata_2018-08-20.csv", check.names = FALSE)

############################################################
# 1. Visits: one row per stool metagenome (1,527 samples)
############################################################

meta <- read.delim("raw/HMP2_metadata.tsv", check.names = FALSE, row.names = 1)
meta$sample_id <- rownames(meta)
taxa <- read.delim("raw/HMP2_taxonomy.tsv", check.names = FALSE, row.names = 1)
taxa <- as.matrix(taxa[meta$sample_id, ])

# Clinical variables recorded at the same visit (full IBDMDB metadata)
clin <- big[match(meta$sample_id, big[["External ID"]]),
            c("External ID", "Immunosuppressants (e.g. oral corticosteroids)",
              "hbi", "sccai", "fecalcal")]

# Dysbiosis score (Lloyd-Price et al. 2019): median Bray-Curtis
# dissimilarity of a sample to non-IBD reference samples from
# other participants.
bc <- as.matrix(vegan::vegdist(taxa, method = "bray"))
ref <- meta$diagnosis == "nonIBD"
dys_score <- vapply(seq_len(nrow(meta)), function(i) {
  keep <- ref & meta$participant_id != meta$participant_id[i]
  median(bc[i, keep])
}, numeric(1))

visits <- meta %>%
  mutate(
    shannon         = vegan::diversity(taxa, index = "shannon"),
    dysbiosis_score = dys_score,
    immunosuppressants = as.integer(clin[["Immunosuppressants (e.g. oral corticosteroids)"]] == "Yes"),
    hbi             = suppressWarnings(as.numeric(clin$hbi)),
    sccai           = suppressWarnings(as.numeric(clin$sccai)),
    fecal_calprotectin = suppressWarnings(as.numeric(clin$fecalcal))
  ) %>%
  transmute(
    sample_id,
    participant_id,
    site       = site_name,
    week       = week_num,
    diagnosis  = factor(diagnosis, levels = c("nonIBD", "UC", "CD")),
    age,
    sex        = factor(sex),
    race,
    probiotic  = as.integer(!grepl("^No", probiotic)),
    immunosuppressants,
    hbi, sccai, fecal_calprotectin,
    antibiotics = as.integer(antibiotics == "Yes"),
    # Dysbiotic: score above the 90th percentile of non-IBD samples (Lloyd-Price et al. 2019).
    # Applied to every participant, so about 10% of non-IBD samples are dysbiotic by definition.
    # (The MaAsLin 3 metadata codes only dysbiosis_CD / dysbiosis_UC, so non-IBD samples appear as "none" there.)
    dysbiosis   = as.integer(dysbiosis_score > quantile(dysbiosis_score[diagnosis == "nonIBD"], 0.9)),
    dysbiosis_score,
    shannon
  ) %>%
  arrange(participant_id, week) %>%
  group_by(participant_id) %>%
  mutate(
    visit            = row_number(),
    weeks_since_prev = week - lag(week),
    antibiotics_prev = lag(antibiotics),
    dysbiosis_prev   = lag(dysbiosis),
    dysbiosis_score_prev = lag(dysbiosis_score),
    shannon_prev     = lag(shannon)
  ) %>%
  ungroup() %>%
  relocate(visit, .after = week)

############################################################
# 2. Species: relative abundances, present in >= 10% of samples
############################################################

keep_sp <- colMeans(taxa > 0) >= 0.10
species <- data.frame(sample_id = rownames(taxa), taxa[, keep_sp],
                      check.names = FALSE)
species <- species[match(visits$sample_id, species$sample_id), ]

############################################################
# 3. Metabolites: stool metabolomics from the same stool sample
#    (metabolomics ID = metagenome ID, with or without a "_P" suffix)
############################################################

mbx <- read.delim("raw/HMP2_metabolites.tsv", check.names = FALSE)
mbx$sample_id <- visits$sample_id[match(mbx$ID, sub("_P$", "", visits$sample_id))]
mbx <- mbx[!is.na(mbx$sample_id), ]
met <- as.matrix(mbx[, setdiff(names(mbx), c("ID", "sample_id"))])
keep_met <- colMeans(met > 0, na.rm = TRUE) >= 0.50     # detected in >= 50%
metabolites <- data.frame(sample_id = mbx$sample_id, met[, keep_met],
                          check.names = FALSE)

############################################################
# 4. Biopsy: host FUT2 expression (log2 CPM) per biopsy
############################################################

htx_raw <- read.delim("raw/host_tx_counts.pcl.tsv", header = FALSE,
                      check.names = FALSE, stringsAsFactors = FALSE)
ann  <- htx_raw[1:13, ]
rownames(ann) <- ann[[1]]
cnt  <- htx_raw[15:nrow(htx_raw), ]
gene <- cnt[[1]]
cnt  <- apply(cnt[, -1], 2, as.numeric)
lib  <- colSums(cnt, na.rm = TRUE)
fut2 <- log2(cnt[gene == "FUT2", ] / lib * 1e6 + 1)

htx_map <- big[big$data_type == "host_transcriptomics",
               c("External ID", "Participant ID")]
htx_map <- htx_map[!duplicated(htx_map[["External ID"]]), ]

biopsy <- data.frame(
  biopsy_id      = unlist(ann["count", -1]),
  biopsy_location = unlist(ann["biopsy_location", -1]),
  inflamed       = as.integer(unlist(ann["is_inflamed", -1]) == "Yes"),
  week           = as.numeric(unlist(ann["week_num", -1])),
  FUT2_log2cpm   = as.numeric(fut2)
)
biopsy$participant_id <- htx_map[["Participant ID"]][match(biopsy$biopsy_id, htx_map[["External ID"]])]
biopsy <- biopsy %>%
  filter(!is.na(participant_id)) %>%
  relocate(participant_id, .after = biopsy_id) %>%
  arrange(participant_id, week, biopsy_location)

############################################################
# 5. Genotype: FUT2 rs516246 (host exome) and ancestry
#    Mapped through host_genome rows only: some External IDs
#    are reused across data types (for example PSM6XC11 is
#    P6010's blood genome sample and P6017's biopsy).
############################################################

vcf <- readLines("raw/FUT2.vcf")
vcf <- vcf[!startsWith(vcf, "##")]
hdr <- strsplit(vcf[1], "\t")[[1]]
rec <- strsplit(vcf[2], "\t")[[1]]
gt  <- sub(":.*", "", rec[10:length(rec)])

geno <- data.frame(genome_id = hdr[10:length(hdr)],
                   rs516246  = c("0/0" = 0L, "0/1" = 1L, "1/1" = 2L)[gt])
gmap <- big[big$data_type == "host_genome", c("External ID", "Participant ID")]
geno$participant_id <- gmap[["Participant ID"]][match(geno$genome_id, gmap[["External ID"]])]

mds <- read.table("raw/HMP.and1kG.MDS.mds", header = TRUE)
geno$ancestry_PC1 <- mds$C1[match(geno$genome_id, mds$IID)]
geno$ancestry_PC2 <- mds$C2[match(geno$genome_id, mds$IID)]

genotype <- geno %>%
  filter(!is.na(participant_id)) %>%          # drops the NA12878 trio controls
  relocate(participant_id) %>%
  arrange(participant_id)

attr(genotype, "variant") <- list(
  rsid = rec[3], chrom = rec[1], pos_GRCh37 = as.integer(rec[2]),
  ref = rec[4], alt = rec[5], coding = "count of ALT alleles")

############################################################
# 6. Data dictionary
############################################################

dictionary <- read.csv(text = '
table,variable,description
visits,sample_id,Stool metagenome sample ID (HMP2 External ID)
visits,participant_id,HMP2 participant ID
visits,site,Recruitment site
visits,week,Study week of the stool sample
visits,visit,Order of the sample within the participant (1 = first)
visits,diagnosis,"nonIBD, UC (ulcerative colitis), or CD (Crohn disease)"
visits,age,Age at consent (years)
visits,sex,Sex
visits,race,Self-reported race
visits,probiotic,Probiotic use in the last 7 days (1 = yes)
visits,immunosuppressants,"Immunosuppressant use, e.g. oral corticosteroids (1 = yes)"
visits,hbi,Harvey-Bradshaw Index of Crohn disease activity (CD only; missing at many visits)
visits,sccai,Simple Clinical Colitis Activity Index (UC only; missing at many visits)
visits,fecal_calprotectin,Fecal calprotectin (ug/g); measured at a subset of visits
visits,antibiotics,Antibiotic use reported at the visit (1 = yes)
visits,dysbiosis,"Dysbiotic sample (1 = yes): dysbiosis_score above the 90th percentile of non-IBD samples (Lloyd-Price et al. 2019); about 10% of non-IBD samples by definition"
visits,dysbiosis_score,Median Bray-Curtis dissimilarity to non-IBD samples from other participants (higher = more dysbiotic)
visits,shannon,Shannon diversity of species relative abundances
visits,weeks_since_prev,Weeks since the participant previous stool sample
visits,antibiotics_prev,antibiotics at the previous sample
visits,dysbiosis_prev,dysbiosis at the previous sample
visits,dysbiosis_score_prev,dysbiosis_score at the previous sample
visits,shannon_prev,shannon at the previous sample
species,sample_id,Matches visits$sample_id
species,(species columns),"Relative abundance (proportion, MetaPhlAn) of species present in at least 10% of samples"
metabolites,sample_id,Matches visits$sample_id (metabolomics from the same stool sample)
metabolites,(metabolite columns),Untargeted LC-MS intensity of named metabolites detected in at least 50% of samples (not normalized)
biopsy,biopsy_id,Biopsy RNA-seq sample ID
biopsy,participant_id,HMP2 participant ID
biopsy,biopsy_location,Biopsy site (for example Rectum or Ileum)
biopsy,inflamed,Biopsy taken from inflamed tissue (1 = yes)
biopsy,week,Study week of the biopsy (most at the baseline colonoscopy)
biopsy,FUT2_log2cpm,FUT2 expression as log2(counts per million + 1)
genotype,participant_id,HMP2 participant ID
genotype,genome_id,Host exome sample ID
genotype,rs516246,"Number of ALT (T) alleles at FUT2 rs516246 (chr19:49206172, GRCh37; C>T): 0, 1, 2"
genotype,ancestry_PC1,First ancestry component (MDS with 1000 Genomes)
genotype,ancestry_PC2,Second ancestry component (MDS with 1000 Genomes)
', strip.white = TRUE)

############################################################
# 7. Save
############################################################

iHMP_causal <- list(visits = as.data.frame(visits), species = species,
                    metabolites = metabolites, biopsy = biopsy,
                    genotype = genotype, dictionary = dictionary)
attr(iHMP_causal, "source") <- paste(
  "HMP2/IBDMDB (Lloyd-Price et al. 2019, Nature 569:655-662).",
  "Built for the Physalia course Causal AI Methods for Computational Biology, 2026.")
saveRDS(iHMP_causal, "out/iHMP_causal.rds")

for (nm in names(iHMP_causal))
  write.csv(iHMP_causal[[nm]], file.path("out", paste0(nm, ".csv")), row.names = FALSE)

cat("visits:     ", nrow(visits), "samples from", n_distinct(visits$participant_id), "participants\n")
cat("species:    ", ncol(species) - 1, "species\n")
cat("metabolites:", nrow(metabolites), "samples x", ncol(metabolites) - 1, "metabolites\n")
cat("biopsy:     ", nrow(biopsy), "biopsies from", n_distinct(biopsy$participant_id), "participants\n")
cat("genotype:   ", nrow(genotype), "participants\n")
