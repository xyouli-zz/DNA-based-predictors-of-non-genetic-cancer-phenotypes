# This script is for calculating gene signature scores based on RNA and segment scores based on gene-level DNA CNA

# Given a gene expression data matrix (gene X sample): edata
# rows are genes in Entrez ID and columns are samples
# run calc_signatures
signature_score <- calc_signatures(edata,"~/data/gene_signatures_20170111.gmt",method = "median")

# find NA signatures
NAsig <- rownames(signature_score[is.na(signature_score[,1]),])
# remove NA signatures
i <- match(NAsig,rownames(signature_score))
signature_score <- signature_score[-i,]

# CD103_Ratio
CD103_pos <- signature_score[rownames(signature_score)=="CD103_Positive_Median_Cancer.Cell.2014_PMID.25446897"]
CD103_neg <- signature_score[rownames(signature_score)=="CD103_Negative_Median_Cancer.Cell.2014_PMID.25446897"]
CD103_ratio <- CD103_pos - CD103_neg # log2 scale division 
signature_score <- rbind(signature_score,CD103_ratio)
rownames(signature_score)[nrow(signature_score)] <- "CD103_Ratio_Cancer.Cell.2014_PMID.25446897"

# differentiation score
diff_centroid <- read.table("~/data/special_gene_signature_training_sets/UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035",sep = '\t',header = T,row.names = 1,check.names = F)
diff_score <- assignDiffScore.dwd(diff_centroid,edata)
signature_score <- rbind(signature_score,diff_score)
rownames(signature_score)[nrow(signature_score)] <- "UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035"

# oncotype DX score
oncotype <- GHI_RS(edata)
signature_score <- rbind(signature_score,diff_score)
rownames(signature_score)[nrow(signature_score)] <- "GHI_RS_Model_NJEM.2004_PMID.15591335"

save(signature_score,file = 'signature_score.rda')

# For special signatures calculated as correlation to predetermined gene centroids
# All traing sets files are included in the ~/data/special_gene_signature_training_sets folder
# For calculation of such special signature, DWD was used to merge current edata with 
# training set, then DWD prediction tool was used to compute pearson/spearman correlation/euclidean distance
# for each sample

# Given a gene-level CNA score matrix (gene X sample): CNdata
# rows are genes in Entrez ID and columns are samples
segment_score <- calc_segments(CNdata,'CNA_segments.gmt',method = 'mean')





