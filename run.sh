#! /bin/bash -e

## It requires out/admixture files and out/plink/DataB (final dataset
## without related and admixed samples).

## Generate main figures
Rscript --vanilla bin/plot-fig1.R
Rscript --vanilla bin/plot-fig2.R
Rscript --vanilla bin/plot-fig3.R
Rscript --vanilla bin/plot-fig4.R
## Generate supplementary material
Rscript --vanilla bin/plot-sup-fig_admixture.R
Rscript --vanilla bin/plot-sup-fig_treemix.R
Rscript --vanilla bin/plot-sup-fig_PCA.R
Rscript --vanilla bin/plot-sup-fig_ROH.R
Rscript --vanilla bin/plot-sup-fig_qpGraph.R