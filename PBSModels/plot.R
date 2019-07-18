library(ggplot2)

# for SNP
#___________________________________________________
snp_df <- read.table(file="ANN_SNP.INFO", header=TRUE, sep="\t",na.strings = "?")

pdf(file="SNP.ANN.distribution.pdf")
qplot(x=AN, data = snp_df, geom = 'density')
qplot(x=DP, data = snp_df, geom = 'density')
qplot(x=QD, data = snp_df, geom = 'density')
qplot(x=MQ, data = snp_df, geom = 'density')
qplot(x=FS, data = snp_df, geom = 'density')
qplot(x=SOR, data = snp_df, geom = 'density')
qplot(x=MQRankSum, data = snp_df, geom = 'density')
qplot(x=ReadPosRankSum, data = snp_df, geom = 'density')
dev.off()

# for INDEL
#___________________________________________________
indel_df <- read.table(file="ANN_INDEL.INFO", header=TRUE, sep="\t",na.strings = "?")

pdf(file="INDEL.ANN.distribution.pdf")
qplot(x=AN, data = indel_df, geom = 'density')
qplot(x=DP, data = indel_df, geom = 'density')
qplot(x=QD, data = indel_df, geom = 'density')
qplot(x=FS, data = indel_df, geom = 'density')
qplot(x=SOR, data = indel_df, geom = 'density')
qplot(x=ReadPosRankSum, data = indel_df, geom = 'density')
dev.off()






