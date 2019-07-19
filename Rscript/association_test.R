# Test for associations between gene signatures and DNA CNA
# Given a signature_score matrix and 
# gene-level CNA matrix CN_score: continuous CNA score for each gene, 
# CN_gain: binary matrix where 1 is copy number gain for a gene and 0 is no gain
# CN_loss: binary matrix where 1 is copy number loss for a gene and 0 is no loss
# Fisher's exact test and spearman correlation test/linear model control for subtypes
# choose a gene signature to test, for example RB-LOH signature

load('BRCA_association_test_data.rda')
pheno <- "UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450"
score <- unlist(signature_score[pheno,])
p_value <- sigCNTest(score, CN_score, CN_gain, CN_loss)

# Benjamini-Hochberg correct p values
p_value_adj <- apply(p_value,2,function(x){return(p.adjust(x,method = "BH"))})
# log10 transform adjusted p values
log_p <- -log(p_value_adj,10)

# plot association landscape
v <- vertical_lines
# for unadjusted associations
p <- log_p[,1:4]
# for subtype-adjusted associations
p <- log_p[,c(5:6,3:4)]

P_Plot(p,main = 'UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450',y1 = 35,y2=50,label_up = seq(0,35,by=10),label_down = seq(10,50,by=10))
