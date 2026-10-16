# Beta-diversity references

Accessed 2026-07-11. Package versions are the versions validated locally; `renv.lock` and the run-level software table are authoritative for a given execution.

## Methods

- Bray JR, Curtis JT. An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs*. 1957;27:325-349. DOI: [10.2307/1942268](https://doi.org/10.2307/1942268). Supports Bray-Curtis dissimilarity.
- Jaccard P. Distribution de la flore alpine dans le bassin des Dranses et dans quelques régions voisines. *Bulletin de la Société Vaudoise des Sciences Naturelles*. 1901;37:241-272. Supports the binary Jaccard coefficient; implementation is verified through `vegan::vegdist`.
- Ružička M. Anwendung mathematisch-statistischer Methoden in der Geobotanik. *Biológia*. 1958;13:647-661. Supports quantitative Jaccard/Ružička. `vegan` quantitative Jaccard is algebraically `2B/(1+B)`, equivalent to Ružička dissimilarity for non-negative profiles.
- Legendre P, Gallagher ED. Ecologically meaningful transformations for ordination of species data. *Oecologia*. 2001;129:271-280. DOI: [10.1007/s004420100716](https://doi.org/10.1007/s004420100716). Supports Hellinger transformation before Euclidean ordination.
- Aitchison J, Barceló-Vidal C, Martín-Fernández JA, Pawlowsky-Glahn V. Logratio analysis and compositional distance. *Mathematical Geology*. 2000;32:271-275. DOI: [10.1023/A:1007529726302](https://doi.org/10.1023/A:1007529726302). Supports CLR and Aitchison distance.
- Martín-Fernández JA, Barceló-Vidal C, Pawlowsky-Glahn V. Dealing with zeros and missing values in compositional data sets using nonparametric imputation. *Mathematical Geology*. 2003;35:253-278. DOI: [10.1023/A:1023866030544](https://doi.org/10.1023/A:1023866030544). Supports principled multiplicative zero replacement.
- Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology*. 2017;8:2224. DOI: [10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224). Supports compositional treatment of microbiome profiles.
- Lozupone C, Knight R. UniFrac: a new phylogenetic method for comparing microbial communities. *Applied and Environmental Microbiology*. 2005;71:8228-8235. DOI: [10.1128/AEM.71.12.8228-8235.2005](https://doi.org/10.1128/AEM.71.12.8228-8235.2005). Supports unweighted UniFrac.
- Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. *Applied and Environmental Microbiology*. 2007;73:1576-1585. DOI: [10.1128/AEM.01996-06](https://doi.org/10.1128/AEM.01996-06). Supports weighted UniFrac.
- Chen J, Bittinger K, Charlson ES, et al. Associating microbiome composition with environmental covariates using generalized UniFrac distances. *Bioinformatics*. 2012;28:2106-2113. DOI: [10.1093/bioinformatics/bts342](https://doi.org/10.1093/bioinformatics/bts342). Supports generalized UniFrac and alpha 0.5 as a balance between endpoints.
- Gower JC. Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika*. 1966;53:325-338. DOI: [10.1093/biomet/53.3-4.325](https://doi.org/10.1093/biomet/53.3-4.325). Supports principal coordinates analysis.
- Legendre P, Anderson MJ. Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. *Ecological Monographs*. 1999;69:1-24. DOI: [10.1890/0012-9615(1999)069[0001:DBRATM]2.0.CO;2](https://doi.org/10.1890/0012-9615%281999%29069%5B0001%3ADBRATM%5D2.0.CO%3B2). Supports Lingoes/Cailliez correction and optional dbRDA.
- Kruskal JB. Nonmetric multidimensional scaling: a numerical method. *Psychometrika*. 1964;29:115-129. DOI: [10.1007/BF02289694](https://doi.org/10.1007/BF02289694). Supports NMDS and stress reporting.
- Anderson MJ. A new method for non-parametric multivariate analysis of variance. *Austral Ecology*. 2001;26:32-46. DOI: [10.1111/j.1442-9993.2001.01070.pp.x](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x). Supports PERMANOVA.
- Anderson MJ. Distance-based tests for homogeneity of multivariate dispersions. *Biometrics*. 2006;62:245-253. DOI: [10.1111/j.1541-0420.2005.00440.x](https://doi.org/10.1111/j.1541-0420.2005.00440.x). Supports `betadisper`/PERMDISP.
- Zhu Z, Satten GA, Mitchell C, Hu YJ. Constraining PERMANOVA and LDM to within-set comparisons by projection improves efficiency of analyses of matched sets of microbiome data. *Microbiome*. 2021;9:133. DOI: [10.1186/s40168-021-01034-9](https://doi.org/10.1186/s40168-021-01034-9). Supports hypothesis-aware within-subject restrictions for longitudinal/matched microbiome analysis.

## Software

- `vegan` 2.7-3: [official reference manual](https://vegandevs.github.io/vegan/reference/). Used for Bray-Curtis, Jaccard/Ružička, PCoA, NMDS, PERMANOVA, and dispersion.
- `permute` 0.9-10: [CRAN package page](https://cran.r-project.org/package=permute). Used to define subject blocks.
- `ape` 5.8-1: [CRAN package page](https://cran.r-project.org/package=ape). Used to read, validate, and prune trees.
- `phyloseq` 1.46.0: [Bioconductor package page](https://bioconductor.org/packages/phyloseq). Used for unweighted and weighted UniFrac when a tree exists.
- `GUniFrac`: [CRAN package page](https://cran.r-project.org/package=GUniFrac). Used for generalized UniFrac when a tree exists.
- `ggplot2` 4.0.3: [official documentation](https://ggplot2.tidyverse.org/). Used for figures.
- `lme4` 2.0-1: [CRAN package page](https://cran.r-project.org/package=lme4). Used for exploratory distance-to-baseline mixed models.

No pairwise-PERMANOVA convenience package is required. Contrasts are implemented internally with the same restricted permutation design and tested against toy and integration cases.
