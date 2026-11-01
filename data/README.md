# iHMP_causal: course dataset

Course data for **Causal AI Methods for Computational Biology** (Physalia, 2 to 6 November 2026). Download `iHMP_causal.rds` into your course folder.

## 1. iHMP_causal.rds (Days 3 and 5)

This dataset comes from the Integrative Human Microbiome Project inflammatory bowel disease cohort (HMP2, IBDMDB): people with Crohn disease, ulcerative colitis, and non-IBD controls, followed for up to a year with stool, biopsy, and blood samples.

## Load it

```r
ihmp <- readRDS("iHMP_causal.rds")
names(ihmp)
visits <- ihmp$visits
```

CSV copies of every table are included for viewing in a spreadsheet.

## Tables

| Table | Rows | Content |
|---|---|---|
| `visits` | 1,527 stool samples from 120 participants | Clinical data at each visit, antibiotic use, dysbiosis, Shannon diversity, and values at the previous visit |
| `species` | 1,527 | Relative abundance of 109 gut species (same rows as `visits`) |
| `metabolites` | 443 | Stool metabolomics, 523 named metabolites, from the same stool samples |
| `biopsy` | 252 biopsies from 90 participants | *FUT2* gene expression in intestinal biopsies |
| `genotype` | 92 participants | *FUT2* variant rs516246 and two ancestry components |
| `dictionary` | | Description of every variable |

Tables link through `sample_id` (stool samples) and `participant_id` (people).

## Source

Lloyd-Price J, Arze C, Ananthakrishnan AN, *et al.* (2019). Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* 569: 655-662. Data portal: <https://ibdmdb.org>.

Species profiles and visit metadata follow the HMP2 files distributed with MaAsLin 3 (<https://github.com/biobakery/maaslin3>). The script `build_iHMP_causal.R` reproduces this dataset from the original files.

## Other course data

Day 5 (TMLE) also uses `RLab5_TMLE.csv`, a small teaching dataset (249 rows: covariates W1 to W4, exposure A, outcome Y) from the Weill Cornell causal inference course.


Day 2 (matching) uses the job training data `lalonde`, which ships with the MatchIt package: `data("lalonde", package = "MatchIt")`.

## Autism microbiota transfer therapy data (Day 4)

`ASD_MTT.rds`: microbiota transfer therapy (MTT) in 18 children with autism (open label; baseline, week 10, 2-year follow-up), 20 typically developing children sampled once, and the stool donors (Kang et al. 2017, *Microbiome*; Kang et al. 2019, *Scientific Reports*). A list with `species`, `fecal_metabolites`, `plasma_metabolites`, `kos` (each with `features`, `metadata`, `features_metadata`) and `symptoms` (CARS, ABC, SRS, GSRS, Vineland, PGI by time point). The same data are in the Zentangler repository (https://github.com/himelmallick/Zentangler).
