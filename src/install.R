# Install BiocManager 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
cran_packages <- c(
  "rlang",
  "dplyr",
  "ggplot2",
  "patchwork",
  "stringr",
  "xgboost"
)

# Bioconductor packages
bioc_packages <- c(
  "curatedMetagenomicData",
  # Pasolli E, Schiffer L, Manghi P, Renson A, Obenchain V, Truong D, Beghini F, 
  # Malik F, Ramos M, Dowd J, Huttenhower C, Morgan M, Segata N, Waldron L (2017). 
  # “Accessible, curated metagenomic data through ExperimentHub.” Nat. Methods, 
  # 14(11), 1023–1024. ISSN 1548-7091, 1548-7105. doi:10.1038/nmeth.4468.
  
  "SummarizedExperiment",
  # Morgan M, Obenchain V, Hester J, Pagès H (2026). SummarizedExperiment: A container 
  # (S4 class) for matrix-like assays. doi:10.18129/B9.bioc.SummarizedExperiment. 
  # R package version 1.42.0, https://bioconductor.org/packages/SummarizedExperiment.
  
  "lefser",
  # Khleborodova A, Gamboa-Tuz S, Ramos M, Segata N, Waldron L, Oh S (2024). 
  # “Lefser: Implementation of metagenomic biomarker discovery tool, LEfSe, in R.” 
  # Bioinformatics, btae707. ISSN 1367-4811. doi:10.1093/bioinformatics/btae707. 
  # https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btae707/7908399.
  
  "CrcBiomeScreen"
  # Li C, Bezbaruah R, Wood H, Gusnanto A (2026). CrcBiomeScreen: An R package for 
  # colorectal cancer screening and microbiome analysis. doi:10.18129/B9.bioc.
  # CrcBiomeScreen. R package version 1.0.0, https://bioconductor.org/packages/CrcBiomeScreen.
)

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Load message
message("All packages installed successfully")
