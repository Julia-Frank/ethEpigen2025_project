# ethEpigen2025_project

# File Overview
All results are included in the main files (project_report.Rmd and project_report.html).  Since preprocessing is quite computationally and memory intensive, we performed it on Euler. All code scripts and resulting data are available in the `preprocessing` folder. The coverage tracks were too large for GitHub, so we uploaded them to [OneDrive](https://1drv.ms/f/c/ad47ab1cb526b6da/EiGJACKkCqhAqehU3_HI-AMBXvwKSJhkNprIwYIos2Fupg?e=tbCstp). If you want to rerun the R script, please download the coverage tracks and place them in `preprocessing`\>`results`\>`tracks`.

# Project proposal

### **Team:** Julia Frank, Sara Leka and Jesslyn Jesslyn

#### 1. What is the topic?

We aim to reproduce the analysis from the study *Application of ATAC-Seq for genome-wide analysis of the chromatin state at single myofiber resolution* (Sahinyan et al., 2022), which used smfATAC-seq to compare chromatin accessibility in uninjured, injured, dystrophic (mdx), and wild-type (WT) mouse myofibers. The mdx mouse is a model for Duchenne Muscular Dystrophy (DMD), a disease marked by progressive muscle degeneration. The study reported distinct accessibility patterns, with mdx associated with muscle structure and WT with metabolic processes. Using computational methods learned in this course, we aim to validate these findings.

#### 2. What data will you be using?

We plan to use the raw sequencing data and reprocess it from scratch. We then aim to compare our resulting normalized genome-wide signal coverage tracks to the provided preprocessed bigWig files to ensure that our preprocessing aligns with theirs.

#### 3. What are the analysis you wish to reproduce or the questions you wish to answer?

Our goal is to compare chromatin accessibility across uninjured, injured, dystrophic (mdx), and wild-type (WT) mouse myofibers. However, if preprocessing proves too time-consuming, we will focus solely on the comparison between mdx and WT myofibers, excluding injured versus uninjured fibers. We aim to identify overlapping and unique accessible regions between the conditions, followed by differential accessibility analysis. Additionally, we plan to perform gene set enrichment and motif enrichment analyses.

#### Links:

The paper: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8901173/#abstract1>\
The data: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173676>
