library(ggplot2)
library(dplyr)
library("stringr", lib.loc="~/miniconda3/envs/Renv/lib/R/library")


# SNP rs1067
pEUR = 0.19980000
pAFR = 0.08930000
nEUR = 1006
nAFR = 1322
# Fischer test
contingency = matrix(
  c(
    nEUR * c(pEUR, 1-pEUR),
    nAFR * c(pAFR, 1-pAFR)
  ),
  byrow = TRUE,
  nrow=2
)
fisher.test(contingency)
# Prop test
prop.test(pEUR*nEUR, nEUR)

#---
# SNP rs3091244
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
snpmart = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
listAttributes(snpmart) %>% filter(str_detect(description, "Chr"))

snps = getBM(
  attributes = c('refsnp_id', 'minor_allele_freq'),
  filters = c('chr_name', 'start', 'end'),
  values = list(8, 100000, 1000000),
  mart = snpmart
)

snps = snps %>% filter(! is.na(minor_allele_freq) & minor_allele_freq > 0.05)
hist(snps$minor_allele_freq)
qplot(snps$minor_allele_freq, geom = "histogram")


# PLINK
# set seed so we all have the same random martrix
set.seed(1)
# create random matrix
size=20
m = 10
genotypes=matrix(rbinom(size*m, size=2, prob=.5), nrow=size)
genotypes
# generate missing values
missing = matrix(sample(c(0, NA), size=size*m, prob=c(.9, .1), replace=T), ncol=m)
missing
# set missing values to genotypes
x = genotypes+missing
qcCol = function(gen, frac){
  gen[, colSums(is.na(gen))/nrow(gen) <= frac]
}
qcRow = function(gen, frac){
  gen[rowSums(is.na(gen))/ncol(gen) <= frac, ]
}
y = qcRow(qcCol(x, .1), .1)
y

# Quality control II
freq = read.table("~/data_polymorphism/hapmap1/plink.frq", header=T)
freq
ggplot(tibble(maf=freq$MAF), aes(maf, ..density..)) +
  geom_histogram(bins=25)


# SNPSTATS
library(snpstats)
assoc = read.table("data_polymorphism/assoc/plink.assoc.logistic", header = T)
assoc
fam = "data_polymorphism/assoc/mudym.fam"
bim = "data_polymorphism/assoc/mudym.bim"
bed = "data_polymorphism/assoc/mudym.bed"
dat = read.plink(bed, bim, fam)
associations = snp.rhs.estimates(affected ~ sex,
                                 family="binomial",
                                 data=dat$fam,
                                 snp.data=dat$genotypes)



