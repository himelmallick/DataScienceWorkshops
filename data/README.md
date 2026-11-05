# Course data

Course data for **Causal AI Methods for Computational Biology** (Physalia, 2 to 6 November 2026).

## Analysis-ready course files

The labs read two files. Download them into your course folder (`install_Causal.R` does this for you).

| File | Days | Content |
|---|---|---|
| `ihmp.rds` | 3 and 5 | One data frame, one row per stool visit (1,527 visits, 120 participants): clinical variables, antibiotics, dysbiosis, Shannon diversity, values at the previous and next visit, *FUT2* genotype and rectal expression, the 109 species at the previous visit, and the difference-in-differences design columns |
| `asd_mtt_mae.rds` | 4 | One `MultiAssayExperiment` with 56 samples (18 autistic children before and after microbiota transfer therapy, 20 typically developing children): symptom scores in the sample data, and species (CLR, relative abundance, counts), species connectivity, fecal and plasma metabolites as aligned layers |

```r
ihmp <- readRDS("ihmp.rds")
mae  <- readRDS("asd_mtt_mae.rds")   # needs the MultiAssayExperiment package
```

Both are built from the source data below.

## Source data

### iHMP_causal.rds

This dataset comes from the Integrative Human Microbiome Project inflammatory bowel disease cohort (HMP2, IBDMDB): people with Crohn disease, ulcerative colitis, and non-IBD controls, followed for up to a year with stool, biopsy, and blood samples.

### Load it

```r
ihmp <- readRDS("iHMP_causal.rds")
names(ihmp)
visits <- ihmp$visits
```

CSV copies of every table are included for viewing in a spreadsheet.

### Tables

| Table | Rows | Content |
|---|---|---|
| `visits` | 1,527 stool samples from 120 participants | Clinical data at each visit, antibiotic use, dysbiosis, Shannon diversity, and values at the previous visit |
| `species` | 1,527 | Relative abundance of 109 gut species (same rows as `visits`) |
| `metabolites` | 443 | Stool metabolomics, 523 named metabolites, from the same stool samples |
| `biopsy` | 252 biopsies from 90 participants | *FUT2* gene expression in intestinal biopsies |
| `genotype` | 92 participants | *FUT2* variant rs516246 and two ancestry components |
| `dictionary` | | Description of every variable |

Tables link through `sample_id` (stool samples) and `participant_id` (people).

### Source

Lloyd-Price J, Arze C, Ananthakrishnan AN, *et al.* (2019). Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* 569: 655-662. Data portal: <https://ibdmdb.org>.

Species profiles and visit metadata follow the HMP2 files distributed with MaAsLin 3 (<https://github.com/biobakery/maaslin3>). The script `build_iHMP_causal.R` reproduces this dataset from the original files.

## Other data

`RLab5_TMLE.csv` is a small teaching dataset (249 rows: covariates W1 to W4, exposure A, outcome Y) from the Weill Cornell causal inference course. The Physalia labs do not use it.


Day 2 (matching) uses the job training data `lalonde`, which ships with the MatchIt package: `data("lalonde", package = "MatchIt")`.

### Autism microbiota transfer therapy data (ASD_MTT.rds)

`ASD_MTT.rds`: microbiota transfer therapy (MTT) in 18 children with autism (open label; baseline, week 10, 2-year follow-up), 20 typically developing children sampled once, and the stool donors (Kang et al. 2017, *Microbiome*; Kang et al. 2019, *Scientific Reports*). A list with `species`, `fecal_metabolites`, `plasma_metabolites`, `kos` (each with `features`, `metadata`, `features_metadata`) and `symptoms` (CARS, ABC, SRS, GSRS, Vineland, PGI by time point). The same data are in the Zentangler repository (https://github.com/himelmallick/Zentangler).
