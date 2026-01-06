#script for running prs
#load libraries

library(tidyverse)

### make lists of txt files (of tibbles) (read in, assign to list)

c1=vroom::vroom("/proj/sens2017538/nobackup/marina/data/clusters_PGS_TS1/clust1FgGt_TS1.txt")
c2=vroom::vroom("/proj/sens2017538/nobackup/marina/data/clusters_PGS_TS1/clust2FgGt_TS1.txt")
c3=vroom::vroom("/proj/sens2017538/nobackup/marina/data/clusters_PGS_TS1/clust3FgGt_TS1.txt")
c4=vroom::vroom("/proj/sens2017538/nobackup/marina/data/clusters_PGS_TS1/clust0FgGt_TS1.txt")

cluster_list = list(clust1FgGt_TS1 = c1,
                    clust2FgGt_TS1 = c2,
                    clust3FgGt_TS1 = c3,
                    clust0FgGt_TS1 = c4)

clusters = names(cluster_list)

#create function for cmd

plink_rsc = "plink --bfile /proj/sens2017538/nobackup/pascal/GIANT_FINNGEN_ANALYSES/genetic_data/gt_fn_SNPs " # can use same files as Pascal
score = "--score /proj/sens2017538/nobackup/marina/data/clusters_PGS_TS1/" #folder where I have SNP files for cluster

cols = "1 2 3 header sum "
out = "--out /proj/sens2017538/nobackup/marina/results/bmipgs_TS1/"

#function
prs_func = function(dat){
  cmd = paste0(plink_rsc,
               score,paste0(dat, ".txt "),cols,
               out,paste0(dat))
  
  system(cmd)
}

#use function to loop over all PRSs

lapply(clusters, function(x) prs_func(x))


## for when running in server
q()
