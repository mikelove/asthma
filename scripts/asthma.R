library(tximport)
library(readr)

coldata <- read_delim("../data/SraRunTable.txt", delim="\t")
files <- file.path("../data/quant",coldata$Run_s,"quant.sf.gz")
names(files) <- coldata$Run_s

if (FALSE) {
  # ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
  library(GenomicFeatures)
  #txdb <- makeTxDbFromGFF("gencode.v26.annotation.gtf.gz")
  #saveDb(txdb, file="gencode.v26.sqlite")
  txdb <- loadDb("gencode.v26.sqlite")

  columns(txdb)
  k <- keys(txdb, "GENEID")
  res <- select(txdb, k, "TXNAME", "GENEID")
  tx2gene <- res[,2:1]
}

load("../data/tx2gene.rda")

# does it work?
# txi <- tximport(files[1], type="salmon", tx2gene=tx2gene)

# ok, now do it for all files
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

library(DESeq2)

# match with GEO title
geo <- read_delim("../data/GEO_table.txt", delim="\t", col_names=FALSE)
coldata$title <- geo$X2[match(coldata$Sample_Name_s, geo$X1)]
coldata$condition <- factor(coldata$disease_state_s)
coldata$treatment <- factor(coldata$treatment_s)

dds <- DESeqDataSetFromTximport(txi, coldata,
                                ~condition + treatment + condition:treatment)

library(magrittr)
levels(dds$condition) <- c("asth","non")
dds$condition %<>% relevel("non")
dds$treatment %<>% relevel("Vehicle")

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, c("treatment","condition"))

dds$id <- substr(dds$title, 1, 3)
id.lvls <- c(dds$id[dds$condition == "non" & dds$treatment == "Vehicle"],
             dds$id[dds$condition == "asth" & dds$treatment == "Vehicle"])

dds$id %<>% factor(id.lvls)
o <- order(dds$condition, dds$treatment, dds$id)
dds <- dds[,o]

as.data.frame(colData(dds)[c("condition","treatment","id")])
all(dds$id == c(rep(id.lvls[1:6], 2),
                rep(id.lvls[7:12], 2)))

dds$id.nested <- factor(rep(1:6,4))

as.data.frame(colData(dds)[c("condition","treatment","id","id.nested")])

design(dds) <- ~condition + condition:id.nested +
  treatment + condition:treatment

rownames(dds) <- make.unique(substr(rownames(dds),1,15))

dds <- dds[rowSums(counts(dds)) > 0,]
keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
res.sort <- res[order(res$pvalue),]

plotMA(res, ylim=c(-5,5))
summary(res)

top.gene <- rownames(res.sort)[1]
plotCounts(dds, top.gene, c("condition","treatment"), transform=FALSE)

library(Homo.sapiens)
Homo.sapiens %>% mapIds(top.gene, "SYMBOL", "ENSEMBL")
go.tab <- Homo.sapiens %>% select(top.gene, "GOID", "ENSEMBL") %>% subset(ONTOLOGY == "BP")
# Homo.sapiens %>% select(go.tab$GOID, "TERM", "GOID")

target <- c("CCL5","CXCL10","CX3CL1","ACKR4","CDHR3")
target.map <- mapIds(Homo.sapiens, target, "ENSEMBL", "SYMBOL")
target.map
match(target.map, rownames(res.sort))
plotCounts(dds, target.map[2], c("condition","treatment"))
plotCounts(dds, target.map[2], c("condition","treatment"), transform=FALSE)

dat <- plotCounts(dds, target.map[2], c("condition","treatment","id.nested"),
                  returnData=TRUE)
library(ggplot2)
ggplot(dat, aes(x=treatment, y=count, col=id.nested, group=id.nested)) +
  geom_point() + geom_smooth(method="lm") +
  scale_y_log10() + 
  facet_wrap(~condition)

#

plotDispEsts(dds)
dispersionFunction(dds)
dmr <- function(x) (.05 + 2.06 / x) * exp(rnorm(length(x),0,.28))
baseMean <- mcols(dds)$baseMean
plot(baseMean, dmr(baseMean), log="xy", ylim=c(1e-8,10))
mean(log2(baseMean))
sd(log2(baseMean))
hist(2^rnorm(1000,mean=4.80,sd=1.66))
set.seed(1)
sim <- makeExampleDESeqDataSet(n=10000, m=12,
                               betaSD=0.5,
                               interceptMean=4.80,
                               interceptSD=1.66,
                               dispMeanRel=dmr)
keep <- rowSums(counts(sim) >= 10) >= 3
table(keep)
sim <- sim[keep,]
sim <- DESeq(sim)
sim.res <- results(sim, independentFiltering=FALSE, cooksCutoff=FALSE)
par(mfrow=c(1,2))
plotMA(res, xlim=c(1,1e6), ylim=c(-5,5))
plotMA(sim.res, xlim=c(1,1e6), ylim=c(-5,5))
sim.dat <- data.frame(sig=sim.res$padj < .1,
                      logmean=cut(mcols(sim)$trueIntercept,c(-1,3.3,6.6,10,20)),
                      abs.lfc=cut(abs(mcols(sim)$trueBeta),c(0,.25,.5,1,3)))
library(dplyr) 
sim.tab <- sim.dat %>% group_by(logmean, abs.lfc) %>% summarize(power=mean(sig))

library(ggplot2)
ggplot(sim.tab, aes(x=abs.lfc, y=power, col=logmean, group=logmean)) + geom_line()
# 

par(mfrow=c(1,1))
de <- rep(c(FALSE,TRUE),c(8000,2000))
hist(2^rnorm(1000,mean=7,sd=2))
set.seed(1)
sim2 <- makeExampleDESeqDataSet(n=10000, m=12,
                                betaSD=ifelse(de, 0.5, 0),
                                interceptMean=6,
                                interceptSD=3,
                                dispMeanRel=dmr)
sim2 <- DESeq(sim2)
sim2.res <- results(sim2)
thresh <- c(1,5,10,15,20)/100
fdr <- sapply(thresh, function(t) {
  sig <- which(sim2.res$padj < t)
  mean(!de[sig])
})
plot(thresh, fdr, ylim=c(0,.3), type="b", col="blue")
abline(0,1)

#

dds2 <- removeResults(dds)
design(dds2) <- ~condition + treatment + condition:id.nested
dds2 <- DESeq(dds2)
resultsNames(dds2)
res2 <- results(dds2, name=c("treatment_HRV16_vs_Vehicle"))
res2 <- results(dds2, contrast=c("treatment","HRV16","Vehicle"))

summary(res2)
plotMA(res2, ylim=c(-10,10))

res2.sort <- res2[order(res2$log2FoldChange, decreasing=TRUE),]

Homo.sapiens %>% mapIds(rownames(res2.sort)[1:40],
                        "SYMBOL", "ENSEMBL")

match(target.map, rownames(res2.sort))

go.tab <- Homo.sapiens %>% select(rownames(res2.sort)[1],
                                  "GO", "ENSEMBL") %>% subset(ONTOLOGY == "BP")

library(GO.db)
go.tab2 <- GO.db %>% select(go.tab$GO, "TERM", "GOID")
substr(go.tab2$TERM, 1, 60)

getTerms <- function(n) {
  go.tab <- Homo.sapiens %>% select(rownames(res2.sort)[n],
                                    "GO", "ENSEMBL") %>% subset(ONTOLOGY == "BP")
  go.tab2 <- Homo.sapiens %>% select(go.tab$GO, "TERM", "GOID")
  substr(go.tab2$TERM, 1, 60)
}

getTerms(2)
getTerms(3)
getTerms(4)
getTerms(5)
getTerms(6)
getTerms(7)
getTerms(8)

#

dds3 <- removeResults(dds)
design(dds3) <- ~condition + treatment

dds3 <- DESeq(dds3)
resultsNames(dds3)
res3 <- results(dds3, lfcThreshold=2)
res3shr <- lfcShrink(dds3, coef=3, res=res3)

plotMA(res3, ylim=c(-15,15))
rs <- rowSums(counts(dds)[,dds$treatment == "Vehicle"])
with(res3[rs < 12,], points(baseMean, log2FoldChange, cex=2, col="dodgerblue"))
plotMA(res3shr, ylim=c(-15,15))
