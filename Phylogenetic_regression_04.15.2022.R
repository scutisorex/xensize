#xdata <- read.csv(file="/Users/elliezack/Desktop/Size Paper/Data and Graphs/Xenarthra Size Paper Data", header = T) #reading in data
library(ggplot2)
library(tidyverse)
library(ape)
library(phytools)
library(nlme)
library(geiger)
library(caper)



candy <- c('#46291b', '#801216', '#cc9a69', '#ad9334', '#c4255b', '#ee97a8', '#94958f', '#2d8f8e', '#ddbdae', '#316f5a')
##          bra       chla        cho         cyc       das         myr       pri         tam         tol
#xdata_log <- xdata %>% add_column(bv.tv.log = log(xdata$BV.TV)) %>% add_column(izl.log = log(xdata$CL)) %>% add_column(tb.n.log = log(xdata$TB.N)) %>% add_column(tb.th.log = log(xdata$TB.TH)) %>% add_column(conn.log = log(xdata$Connectivity)) %>% add_column(conn.d.log = log(xdata$Conn.D)) %>% add_column(csa.log = log(xdata$CSA)) %>% add_column(gc.log = log(xdata$GC)) %>% add_column(mil.log = log(xdata$MIL)) %>% add_column(sld.log = log(xdata$SLD))
data <- read.csv("PGLS_Regression_Data.csv", header = T)
#data_log <- data %>% add_column(bv.tv.log = log(data$BV.TV)) %>% add_column(proxy.log = log(data$mass_proxy)) %>% add_column(tb.n.log = log(data$TB.N)) %>% add_column(tb.th.log = log(data$TB.TH)) %>% add_column(conn.log = log(data$Connectivity)) %>% add_column(conn.d.log = log(data$Conn.D)) %>% add_column(csa.log = log(data$CSA)) %>% add_column(gc.log = log(data$GC)) %>% add_column(mil.log = log(data$MIL)) %>% add_column(sld.log = log(data$SLD))%>% add_column(mass.log = log(data$mass)) %>% add_column(cw.log = log(data$CW))%>% add_column(izl.log = log(data$IZL))
#data_log <- data_log %>% arrange(desc(row_number()))

######## Note: this command does the same thing as the above chunk of code, with a lot fewer opportunities for error, and therefore less debugging! It does change some of the column names but I fixed them in the versions for ps1 below.
data_log <- data %>% 
  mutate_at(c(5:17), log) %>% 
  arrange(desc(row_number()))

##### Data Prep ######
#ps1
ps1 <- data_log %>% 
  filter(pos.for.analysis == "ps1") %>% 
  arrange(desc(row_number()))
rownames(ps1) <- ps1$Taxon



#ps2
ps2 <- data_log %>% 
  filter(pos.for.analysis == "ps2")%>% 
  arrange(desc(row_number()))
rownames(ps2) <- ps2$Taxon


#ps3
ps3 <- data_log %>% 
  filter(pos.for.analysis == "ps3")%>% 
  arrange(desc(row_number()))
rownames(ps3) <- ps3$Taxon


#ps4
ps4 <- data_log %>% 
  filter(pos.for.analysis == "ps4")%>% 
  arrange(desc(row_number()))
rownames(ps4) <- ps4$Taxon


#ps5
ps5 <- data_log %>% 
  filter(pos.for.analysis == "ps5") %>% 
  arrange(desc(row_number()))
rownames(ps5) <- ps5$Taxon


#ps6
ps6 <- data_log %>% 
  filter(pos.for.analysis == "ps6") %>% 
  arrange(desc(row_number()))
rownames(ps6) <- ps6$Taxon




#Creating the phylogeny
guys <- read.csv(file = "specieswithecology.csv", header = T) %>% 
  rename(taxon = ï..taxon) # I had to rename the column with the taxa in it for some reason, my computer didn't read it right!

rownames(guys) <- guys$taxon
bot <- read.tree(file = "Gibbsetal2015_tree.newick.txt")
plot(bot)
xen <- bot$tip.label

comp <- name.check(bot, guys)

ins <- drop.tip(bot, comp$tree_not_data)
write.tree(ins, file = 'OG_Tree_species.txt')
plot(ins)
name.check(ins, ps1)


# generalized linear model

# Steph's edits for ps1:

# Run the PGLS: 
#-------------------NOTE!!! You couldn't see the lines because your variables were in the wrong order!-------------------------#

bvtv1<-gls(BV.TV~mass_proxy, data=ps1, correlation=corBrownian(1, phy = ins, form = ~Taxon))# You need the "form" term because of some update to the package in the last couple years. It gives you a warning if you don't have it. And it refers back to the data, so I think you have to include it in the gls command itself in full so that it can refer to it. It works fine inside your loop (see below).
summary(bvtv1)

# Run a linear model without the phy correlation:
nophybvtv <- gls(BV.TV~mass_proxy, data=ps1)

# plot the log-transformed variables:
plot(BV.TV~mass_proxy, data = data_log, pch = 16)

# plot the regular GLS, teal
abline(nophybvtv,lwd = 7, col = "#edb742")

# plot the PGLS, tangerine!
abline(bvtv1, lwd = 7, col = "#0e8c90")



bvtv2<-gls(CL~BV.TV, data=ps2, correlation=bm)
summary(bvtv2)
bvtv3<-gls(CL~BV.TV, data=ps3, correlation=bm)
summary(bvtv3)
bvtv4<-gls(CL~BV.TV, data=ps4, correlation=bm)
summary(bvtv4)
bvtv5<-gls(CL~BV.TV, data=ps5, correlation=bm)
summary(bvtv5)
bvtv6<-gls(CL~BV.TV, data=ps6, correlation=bm)
summary(bvtv6)

# pgls <- function(x){
#   metric <- c('BV.TV', 'TB.N', 'TB.TH', 'GC', 'CSA', 'Conn.D', 'Connectivity')
#   pos <- c('ps1', 'ps2', 'ps3', 'ps4', 'ps5', 'ps6')
#   for (i in 1:6)
#     ps <- x %>% filter(pos.for.analysis == pos[i])
#      for (j in 1:7)
#        y <- ps %>% select(matches(metric[j]))
#        # gls <- gls(CL ~ y, data=ps, correlation = bm)
#        # print(summary(gls))
# }


pos <- c('ps1', 'ps2', 'ps3', 'ps4', 'ps5', 'ps6')
for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(BV.TV~mass_proxy, data = ps, correlation = corBrownian(1, phy = ins, form = ~Taxon)) #note I added the "form" term here and switched the order of variables
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.n.log, data = ps, correlation = bm)
  print(summary(gls))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.th.log, data = ps, correlation = bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~gc.log, data = ps, correlation = bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~csa.log, data = ps, correlation = bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mil.log, data = ps, correlation = bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~conn.d.log, data = ps, correlation = bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~conn.log, data = ps, correlation = bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mass.log, data = ps, correlation = bm)
  print(summary(gls))
}


###### lambda ######
# ins$edge.length
# bv.tv.log1 <- data_log %>% filter(pos.for.analysis == 'ps1') %>% select(bv.tv.log, Taxon)
# phylosig(ins, bv.tv.log1, method = 'lambda', test = TRUE)

## geomorph and K mult for multivariate analyses (look for K mult and blomberg's K)


######## Pagel's Lambda ########

# fitPagel <- gls(matur.L ~ age.mat, correlation=corPagel(value=0.8, phy=tree3), data=dat3)
# intervals(fitPagel, which="var-cov")

pl_bm <- corPagel(0, ins, fixed = TRUE)
pl_phy <- corPagel(1, ins, fixed = TRUE)
pl <- corPagel(0.5, ins, fixed = TRUE)

ps <- data_log %>% filter(pos.for.analysis == pos[1])
gls <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl)
gls_bm <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl_bm)
gls_phy <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl_phy)

anova(gls, gls_phy)
anova(gls, gls_bm)

lambda <- seq(0, 1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                  correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                  data = ps)))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                       lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl)
  gls_bm <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl_bm)
  gls_phy <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl_phy)
  
  anova(gls, gls_phy)
  anova(gls, gls_bm)
  
  lambda <- seq(0, 1, length.out = 500)
  lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                    correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                    data = ps)))
  plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                         lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.n.log, data = ps, correlation = pl)
  gls_bm <- gls(proxy.log~tb.n.log, data = ps, correlation = pl_bm)
  gls_phy <- gls(proxy.log~tb.n.log, data = ps, correlation = pl_phy)
  
  anova(gls, gls_phy)
  anova(gls, gls_bm)
  
  lambda <- seq(0, 1, length.out = 500)
  lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                    correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                    data = ps)))
  plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                         lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.th.log, data = ps, correlation = pl)
  gls_bm <- gls(proxy.log~tb.th.log, data = ps, correlation = pl_bm)
  gls_phy <- gls(proxy.log~tb.th.log, data = ps, correlation = pl_phy)
  
  anova(gls, gls_phy)
  anova(gls, gls_bm)
  
  lambda <- seq(0, 1, length.out = 500)
  lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                    correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                    data = ps)))
  plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                         lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~gc.log, data = ps, correlation = pl)
  gls_bm <- gls(proxy.log~gc.log, data = ps, correlation = pl_bm)
  gls_phy <- gls(proxy.log~gc.log, data = ps, correlation = pl_phy)
  
  anova(gls, gls_phy)
  anova(gls, gls_bm)
  
  lambda <- seq(0, 1, length.out = 500)
  lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                    correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                    data = ps)))
  plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                         lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~csa.log, data = ps, correlation = pl)
  gls_bm <- gls(proxy.log~csa.log, data = ps, correlation = pl_bm)
  gls_phy <- gls(proxy.log~csa.log, data = ps, correlation = pl_phy)
  
  anova(gls, gls_phy)
  anova(gls, gls_bm)
  
  lambda <- seq(0, 1, length.out = 500)
  lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                    correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                    data = ps)))
  plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                         lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mil.log, data = ps, correlation = pl)
  gls_bm <- gls(proxy.log~mil.log, data = ps, correlation = pl_bm)
  gls_phy <- gls(proxy.log~mil.log, data = ps, correlation = pl_phy)
  
  anova(gls, gls_phy)
  anova(gls, gls_bm)
  
  lambda <- seq(0, 1, length.out = 500)
  lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                    correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                    data = ps)))
  plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                         lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~conn.d.log, data = ps, correlation = pl)
  gls_bm <- gls(proxy.log~conn.d.log, data = ps, correlation = pl_bm)
  gls_phy <- gls(proxy.log~conn.d.log, data = ps, correlation = pl_phy)
  
  anova(gls, gls_phy)
  anova(gls, gls_bm)
  
  lambda <- seq(0, 1, length.out = 500)
  lik <- sapply(lambda, function(lambda) logLik(gls(proxy.log ~ bv.tv.log,
                                                    correlation = corPagel(value = lambda, phy = ins, fixed = TRUE),
                                                    data = ps)))
  plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                         lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls_bm <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls_phy <- gls(proxy.log~bv.tv.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.n.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.n.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.th.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.th.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~gc.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~gc.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~csa.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~csa.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mil.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mil.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~conn.d.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~conn.d.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mass.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mass.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}


######## Blomberg's K ########
bk_bm <- corBlomberg(1, ins, fixed = TRUE)
bk_phy <- corBlomberg(5, ins)


pos <- c('ps1', 'ps2', 'ps3', 'ps4', 'ps5', 'ps6')

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls_bm <- gls(proxy.log~bv.tv.log, data = ps, correlation = bk_bm)
  print(summary(gls))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls_phy <- gls(proxy.log~bv.tv.log, data = ps, correlation = bk_phy)
  print(summary(gls))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.n.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.n.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.th.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~tb.th.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}


for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~gc.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~gc.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~csa.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~csa.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mil.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mil.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~conn.d.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~conn.d.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mass.log, data = ps, correlation = pl_bm)
  print(summary(gls))
}

for (i in 1:6){
  ps <- data_log %>% filter(pos.for.analysis == pos[i])
  gls <- gls(proxy.log~mass.log, data = ps, correlation = pl_phy)
  print(summary(gls))
}


