library(consensus)
library(ggplot2)
library(GenomicRanges)
library(effsize)
setwd("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform")
load("Full_complement/Mvalues_FINAL.RData")

stopifnot(all(rownames(EPICv1_Ms) %in% rownames(WGBS_Ms)))
singlepos <- rownames(EPICv2_Ms)[new_manifest$posrep==""]

EPICv1_Ms_nodup <- EPICv1_Ms[singlepos,]
EPICv2_Ms_nodup <- EPICv2_Ms[singlepos,]
WGBS_Ms_nodup <- WGBS_Ms[singlepos,]

df <- data.frame(sample=colnames(EPICv1_Ms_nodup))
df$col <- ifelse(grepl("^FD", df$sample), "blue", ifelse(grepl("SYN", df$sample), "red", "white"))
df$text <- ""
df$text[c(1,2,9,10)] <- c("PrEC", "LNCaP", "MCF7", "TAMR")

#Figure 5a, 5b, 5c
#PCAs
#EPICv1
pca <- prcomp(t(EPICv1_Ms_nodup))
plot(pca$x[,1:2], xlab="", ylab="", pch=16, col=df$col, cex=2)
text(pca$x[,1:2], labels = df$text, cex=1.2)
lam <- pca$sdev
v <- lam^2
v <- v/sum(v)
k1 <- round(100*v[1], 2)
k2 <- round(100*v[2], 2)
title(main=paste0("EPICv1 PCA, ", k1+k2, "% of variance explained"), xlab=paste0(k1, "% of variance explained"), ylab=paste0(k2, "% of variance explained"))
legend("topleft", c("Prostate samples", "PDX samples"), col=c("red", "blue"), pch=16, bty="n")

pca <- prcomp(t(EPICv2_Ms_nodup))
plot(pca$x[,1:2], xlab="", ylab="", pch=16, col=df$col, cex=2)
text(pca$x[,1:2], labels = df$text, cex=1.2)
lam <- pca$sdev
v <- lam^2
v <- v/sum(v)
k1 <- round(100*v[1], 2)
k2 <- round(100*v[2], 2)
title(main=paste0("EPICv2 PCA, ", k1+k2, "% of variance explained"), xlab=paste0(k1, "% of variance explained"), ylab=paste0(k2, "% of variance explained"))
legend("topleft", c("Prostate samples", "PDX samples"), col=c("red", "blue"), pch=16, bty="n")

pca <- prcomp(t(WGBS_Ms_nodup))
pca$x[,2] <- -pca$x[,2]
plot(pca$x[,1:2], xlab="", ylab="", pch=16, col=df$col, cex=2)
text(pca$x[,1:2], labels = df$text, cex=1.2)
lam <- pca$sdev
v <- lam^2
v <- v/sum(v)
k1 <- round(100*v[1], 2)
k2 <- round(100*v[2], 2)
title(main=paste0("WGBS PCA, ", k1+k2, "% of variance explained"), xlab=paste0(k1, "% of variance explained"), ylab=paste0(k2, "% of variance explained"))
legend("topleft", c("Prostate samples", "PDX samples"), col=c("red", "blue"), pch=16, bty="n")



#consensus

pfnames=c("EPICv1", "EPICv2", "WGBS")

EPICmm <- MultiMeasure(names=pfnames, data=list(EPICv1_Ms_nodup, EPICv2_Ms_nodup, WGBS_Ms_nodup))
fit <- fitConsensus(EPICmm)

#Marginals
plotMarginals(fit, "average", pal=brewer.pal(n=3, name="Paired"))
plotMarginals(fit, "sensitivity", pal=brewer.pal(n=3, name="Paired"), xlim=c(-0.5, 3.5))
plotMarginals(fit, "precision", pal=brewer.pal(n=2, name="Paired"))

cohen.d(fit@b_i[,2], fit@b_i[,1], paired = T)

#Cohen's d

#d estimate: 0.1386603 (negligible)
#95 percent confidence interval:
#    lower     upper 
#0.1372061 0.1401146

cohen.d(fit@b_i[,3], fit@b_i[,1], paired = T)

#Cohen's d

#d estimate: 1.569712 (large)
#95 percent confidence interval:
#   lower    upper 
#1.562145 1.577279 

cohen.d(fit@b_i[,3], fit@b_i[,2], paired = T)

#Cohen's d

#d estimate: 1.481058 (large)
#95 percent confidence interval:
#   lower    upper 
#1.473726 1.488390 

cohen.d(log(fit@d_i[,2]), log(fit@d_i[,1]), paired = T)

#Cohen's d

#d estimate: 0.03166032 (negligible)
#95 percent confidence interval:
#     lower      upper 
#0.02949540 0.03382523 

cohen.d(log(fit@d_i[,3]), log(fit@d_i[,1]), paired = T)

#Cohen's d

#d estimate: 1.695421 (large)
#95 percent confidence interval:
#   lower    upper 
#1.693340 1.697501 

cohen.d(log(fit@d_i[,3]), log(fit@d_i[,2]), paired = T)

#Cohen's d

#d estimate: 1.662952 (large)
#95 percent confidence interval:
#   lower    upper 
#1.660935 1.664969 

#Figure S12a
#Average vs sensitivity
par(mfrow=c(1, 3))
smoothScatter(fit@a_i[,1], fit@b_i[,1], xlab="Intercept (average)", ylab="Sensitivity", ylim=c(-1, 2), main="EPICv1")
smoothScatter(fit@a_i[,2], fit@b_i[,2], xlab="Intercept (average)", ylab="Sensitivity", ylim=c(-1, 2), main="EPICv2",
              colramp = colorRampPalette(brewer.pal(9, "Purples")))
smoothScatter(fit@a_i[,3], fit@b_i[,3], xlab="Intercept (average)", ylab="Sensitivity", ylim=c(-0.5, 3.5), main="WGBS",
              colramp = colorRampPalette(brewer.pal(9, "Greens")))

library(RColorBrewer)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), 2)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste0(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, substitute(paste(italic('r'), " = ", r)), cex=2)
}

#Figure S12b
pairs(fit@a_i, panel = function(...) smoothScatter(..., nrpoints = 0, add = TRUE, colramp = colorRampPalette(brewer.pal(9, "Purples"))), gap = 0.2, labels = c("EPICv1", "EPICv2", "WGBS"), 
      main=expression(Joint~italic('a')['i']~-~methylation~data), upper.panel = panel.cor)

#Difference between platforms
df <- data.frame(Sensitivity=as.numeric(fit@b_i), Precision=as.numeric(fit@d_i), Platform=c(rep("EPICv1", nrow(fit@b_i)),
                                                                 rep("EPICv2", nrow(fit@b_i)),
                                                                 rep("WGBS", nrow(fit@b_i))))

#Figure 5d
ggplot(df, aes(x=factor(Platform), y=Sensitivity, fill=Platform)) +
  geom_violin() + 
  stat_summary(mapping = aes(x = Platform, y = Sensitivity, fill=Platform),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median)

#Figure 5e
ggplot(df, aes(x=factor(Platform), y=Precision, fill=Platform)) +
  geom_violin() + scale_y_continuous(trans = "log", labels=function(x) sprintf("%.2f", x)) +
  stat_summary(mapping = aes(x = Platform, y = Precision, fill=Platform),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median)


#Probe Type
probetype <- new_manifest$Infinium_Design_Type
names(probetype) <- rownames(EPICv2_Ms)
probetype <- probetype[singlepos]
df <- data.frame(Probe.type=rep(probetype, 2), Platform=c(rep("EPICv1", 586916), rep("EPICv2", 586916)),
                 Average=c(fit@a_i[,"EPICv1"], fit@a_i[,"EPICv2"]),
                 Sens=c(fit@b_i[,"EPICv1"], fit@b_i[,"EPICv2"]),
                 Prec=c(fit@d_i[,"EPICv1"], fit@d_i[,"EPICv2"]))


#Figure 5f
ggplot(df, aes(x=factor(Probe.type), y=Sens, fill=Platform)) +
  geom_violin() + scale_fill_manual(values=brewer.pal(n=2, name="Paired")) +
  ggtitle("Sensitivity by probe type") +
  labs(x="Infinium Probe type") + 
  stat_summary(mapping = aes(x = Probe.type, y = Sens, fill=Platform),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, position = position_dodge(width=0.9))
vp2

#Figure 5g
vp3 <-  ggplot(df, aes(x=factor(Probe.type), y=Prec, fill=Platform)) +
  geom_violin() + scale_fill_manual(values=brewer.pal(n=2, name="Paired")) +
  ggtitle("Precision by probe type") +
  labs(x="Infinium Probe type") + scale_y_continuous(trans = "log") +
  stat_summary(mapping = aes(x = Probe.type, y = Prec, fill=Platform),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, position = position_dodge(width=0.9))
vp3
#Effect sizes
dfE2 <- df[df$Platform=="EPICv2",]
cohen.d(Average~Probe.type, data=dfE2)
#d estimate: -0.9874221 (large)

cohen.d(Sens~Probe.type, data=dfE2)
#d estimate: -0.3769124 (small)

cohen.d(log(Prec)~Probe.type, data=dfE2)
#d estimate: -0.3576897 (small)

dfE1 <- df[df$Platform=="EPICv1",]
cohen.d(Average~Probe.type, data=dfE1)
#d estimate: -1.031964 (large)

cohen.d(Sens~Probe.type, data=dfE1)
#d estimate: -0.2721001 (small)

cohen.d(log(Prec)~Probe.type, data=dfE1)
#d estimate: -0.2898113 (small)

