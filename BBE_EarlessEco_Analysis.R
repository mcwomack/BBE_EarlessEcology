#load libraries
library(tidyverse)
library(phytools)
library(geiger)
library(cowplot)
library(ggdist)
library(colorspace)
library(caper)
library(grid)
library(viridis)

setwd("~/Documents/Documents - Molly’s MacBook Pro/Pubs/SubmittedOrRevising/Earless_Ecology")

#load eardata 
eardata<-read.csv('eardata.csv')

#choose and load phylogeny 
AnuranTree<-read.tree('Data/TreePL-Rooted_Anura_bestTree.tre')  

#check how many genera lack ear data 
length(unique(eardata$Genus_species))
length(unique(eardata$family)) 
length(unique(eardata$genus)) 
setdiff(allsp$genus,eardata$genus) #identify any species from ear dataset you cannot match


##### FUNCTIONS USED DOWNSTREAM #####
#prepare data for caper PGLS analysis
caperdata <- function(data, tree) { 
  caperdata <-data[tree$tip.label, ] #Align data and phylogeny labels
  geocapertree <-makeNodeLabel(tree, method = "number", prefix = "Node")
  return(comparative.data(phy=geocapertree, data=caperdata, names.col=Genus_species, vcv=TRUE, na.omit=FALSE, warn.dropped=TRUE))
}

#trim phylogenetic tree to species within a dataset
trimtree <- function(data, tree) { 
  name.check(tree,data)->overlap
  #overlap
  drop.tip(tree,overlap$tree_not_data)
  #plot(tree)
}

#Remove excess data from species not found in tree
trimdata<-function(data, tree) {
  name.check(tree,data)->overlap
  data<-data[ ! row.names(data) %in% overlap$data_not_tree, ]
  data[] <- lapply(data, function(x) if(is.factor(x)) factor(x) else x)
  as.data.frame(data)
}

#get proportions for proportion plots
getprops <- function(data, x) {
  xdata<-data %>% group_by(ear,{{x}}) %>%
    summarise(earedcount = sum(ear=="Y"),
              earlesscount = sum(ear=="N"))
  xdata$earlessprop<-xdata$earlesscount/(sum(xdata$earlesscount))
  xdata$earedprop<-xdata$earedcount/(sum(xdata$earedcount))
  xdata$n<-xdata$earedcount+xdata$earlesscount
  xdata$prop<-xdata$earlessprop + xdata$earedprop
  cbind(xdata)
}

#proportion plot
propplot <- function(data, x, y, z, colrpal) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{z}})) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colrpal)+
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}


##### MAKES EARLESS ONLY DATASET FOR DATA EXPLORATION #####

fullearless <- eardata %>% filter(ear == 'N')
fulleared <- eardata %>% filter(ear == 'Y')

row.names(fullearless)<-fullearless$Genus_species
eartree<-trimtree(eardata, AnuranTree)
earlesstree<-trimtree(fullearless, AnuranTree)
treeearless<-trimdata(fullearless, earlesstree)

length(fullearless[,1])
length(treeearless[,1])

#explore percentages of earless species with certain conditions
earlessgeo <- fullearless %>% filter(!is.na(alt_min))
earlessmicro <- fullearless %>% filter(!is.na(Microhabitat))
earlessactiv <- fullearless %>% filter(!is.na(ActiveTime), ActiveTime != '__', ActiveTime != 'D_C_N', ActiveTime != 'D__N', ActiveTime != '_C_')
earlessBS <- fullearless %>% filter(!is.na(Collapsed_SVL), Collapsed_SVL > 0)
earlessdf <- fullearless %>% filter(!is.na(df), df > 0)

write.csv(treeearless,"tree_earless.csv")

################################################## GEOGRAPHY ANALYSIS ###########################

##### RUN GEO PGLS #####
geodata <- eardata %>% filter(!is.na(alt_min))
length(geodata[,1])
table(geodata$ear)

#generate range column
geodata$latitude_range<-geodata$latitude_max-geodata$latitude_min
geodata$alt_range<-geodata$alt_max-geodata$alt_min

row.names(geodata)<-geodata$Genus_species
geotree<-trimtree(geodata, AnuranTree)
geodata<-trimdata(geodata, geotree)
table(geodata$ear)

setdiff(fullearless$Genus_species,AnuranTree$tip.label)

geodata <- geodata %>% mutate(ear = if_else(ear== "Y", 1, 0))
geodata$latitude_range_1 <- geodata$latitude_range + 1
geodata$alt_min_78 <- geodata$alt_min + 78
geocaperdata<-caperdata(geodata, geotree)

altmin.pgls<-pgls(ear~sqrt(alt_min_78), data= geocaperdata, lambda="ML")
altmax.pgls<-pgls(ear~sqrt(alt_max), data= geocaperdata, lambda="ML")
altmean.pgls<-pgls(ear~sqrt(alt_mean), data= geocaperdata, lambda="ML")
altrange.pgls<-pgls(ear~sqrt(alt_range), data= geocaperdata, lambda="ML")
latmin.pgls<-pgls(ear~latitude_min, data= geocaperdata, lambda="ML")
latmax.pgls<-pgls(ear~latitude_max, data= geocaperdata, lambda="ML")
latmean.pgls<-pgls(ear~latitude_mean, data= geocaperdata, lambda="ML")
latrange.pgls<-pgls(ear~log(latitude_range_1), data= geocaperdata, lambda="ML")

par(mfrow=c(2,2))
plot(altmin.pgls) 
plot(altmean.pgls) 
plot(altmax.pgls)
plot(altrange.pgls) 
plot(latmin.pgls)
plot(latmean.pgls)
plot(latmax.pgls) 
plot(latrange.pgls)

summary(altmin.pgls)
summary(altmean.pgls)
summary(altmax.pgls)
summary(altrange.pgls)
summary(latmin.pgls)
summary(latmean.pgls)
summary(latmax.pgls)
summary(latrange.pgls)

pdf(file = "Supp_Fig1_Alt_Min.pdf", width = 6,  height = 6)
par(mfrow=c(2,2))
plot(altmin.pgls)
dev.off()

##### PLOT GEOGRAPHIC DATA #####
geodata <- eardata %>% filter(!is.na(alt_min))
length(geodata[,1])
table(geodata$ear)

#generate range column
geodata$latitude_range<-geodata$latitude_max-geodata$latitude_min
geodata$alt_range<-geodata$alt_max-geodata$alt_min

row.names(geodata)<-geodata$Genus_species
#geotree<-trimtree(geodata, AnuranTree)
#geodata<-trimdata(geodata, geotree)

geoplot1 <- function(data, x, y) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{x}})) +
    geom_violin(alpha = 0.5) +
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
    scale_fill_manual(values=c("blue","orange")) +
    theme_classic()
}
geoplot2 <- function(data, x, y, sumstat) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{x}})) + 
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .25, outlier.shape = NA) +
    geom_point(size = 1, alpha = .2, position = position_jitter(seed = 1, width = .1)) + 
    annotation_custom(grobTree(textGrob(paste(" F = ", round(sumstat$fstatistic[1], digits = 2), "\n","adj. R2 = ", round(sumstat$adj.r.squared, , digits = 4),"\n", "p = ", round(sumstat$coefficients[2,4], digits = 3)), x=0.4,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=12, fontface="italic")))) +
    #coord_cartesian(xlim = c(1.2, NA), clip = "off") +
    theme_classic() + 
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_manual(values = pal, guide = "none") 
}

geodata$ear <- factor(geodata$ear, levels = c('Y', 'N'), ordered = T)
pal <- c("#FF8C00", "#2D7297")

ggalt1<-geoplot2(geodata, ear, alt_min, summary(altmin.pgls))
ggalt3<-geoplot2(geodata, ear, alt_max, summary(altmax.pgls))
ggalt2<-geoplot2(geodata, ear, alt_mean, summary(altmean.pgls))
ggalt4<-geoplot2(geodata, ear, alt_range, summary(altrange.pgls))
gglat1<-geoplot2(geodata, ear, latitude_min, summary(latmin.pgls))
gglat3<-geoplot2(geodata, ear, latitude_max, summary(latmax.pgls))
gglat2<-geoplot2(geodata, ear, latitude_mean, summary(latmean.pgls))
gglat4<-geoplot2(geodata, ear, latitude_range, summary(latrange.pgls))
plot_grid(gglat1, gglat2, gglat3, gglat4, ggalt1, ggalt2, ggalt3, ggalt4, labels = c('A', 'B','C','D','E','F'))

pdf(file = "Fig1_GeoPlots_.pdf", width = 16,  height = 8)
plot_grid(gglat1, gglat2, gglat3, gglat4, ggalt1, ggalt2, ggalt3, ggalt4, labels = c('B','C','D','E','F','G','H','I'), ncol = 4)
dev.off()


##### SUPPLEMENTAL EXAMINATION OF RAPOPORT'S RULE #####
#trim geodata to species below 25 degrees latitude 
rapdata25 <- geodata %>% filter(abs(latitude_max) < 25)

row.names(rapdata25)<-rapdata25$Genus_species
raptree25<-trimtree(rapdata25, AnuranTree)
table(rapdata25$ear)

#prepare dataset for cape PGLS 
rapcaperdata25<-caperdata(rapdata25, raptree25)

rap25_latrange.pgls<-pgls(ear~log(latitude_range_1), data= rapcaperdata25, lambda="ML")
rap25_altrange.pgls<-pgls(ear~sqrt(alt_range), data= rapcaperdata25, lambda="ML")

summary(rap25_latrange.pgls)
summary(rap25_altrange.pgls)

rapdata15 <- geodata %>% filter(abs(latitude_max) < 15)

row.names(rapdata15)<-rapdata15$Genus_species
raptree15<-trimtree(rapdata15, AnuranTree)
table(rapdata15$ear)

#prepare dataset for cape PGLS 
rapcaperdata15<-caperdata(rapdata15, raptree15)

rap15_latrange.pgls<-pgls(ear~log(latitude_range_1), data= rapcaperdata15, lambda="ML")
rap15_altrange.pgls<-pgls(ear~sqrt(alt_range), data= rapcaperdata15, lambda="ML")

summary(rap15_latrange.pgls)
summary(rap15_altrange.pgls)

pal <- c("#2D7297","#FF8C00")

gg1<-geoplot2(rapdata25, as.factor(ear), latitude_range, summary(rap25_latrange.pgls))
gg2<-geoplot2(rapdata15, as.factor(ear), latitude_range, summary(rap15_latrange.pgls))
gg3<-ggplot(geodata, aes(x = abs(latitude_max), y = latitude_range, color = as.factor(ear))) +
  geom_point(size = 3, alpha = .7, position = position_jitter(seed = 1, width = 0.2)) +
  geom_smooth(method = "lm", fill = NA) +
  scale_color_manual(values=pal) +
  theme_classic()
gg4<-geoplot2(rapdata25, as.factor(ear), alt_range, summary(rap25_altrange.pgls))
gg5<-geoplot2(rapdata15, as.factor(ear), alt_range, summary(rap15_altrange.pgls))
gg6<-ggplot(geodata, aes(x = abs(latitude_max), y = alt_range, color = as.factor(ear))) +
  geom_point(size = 3, alpha = .7, position = position_jitter(seed = 1, width = 0.2)) +
  geom_smooth(method = "lm", fill = NA) +
  scale_color_manual(values=pal) +
  theme_classic()


pdf(file = "SuppFig1_RapoportLat_.pdf", width = 8,  height = 8)
toprow <- plot_grid(gg1, gg2, labels = c('A', 'B'), label_size = 12)
plot_grid(toprow, gg3, labels = c('', 'C'), label_size = 12, ncol = 1)
dev.off()

pdf(file = "SuppFig2_RapoportAlt_.pdf", width = 8,  height = 8)
toprow <- plot_grid(gg4, gg5, labels = c('A', 'B'), label_size = 12)
plot_grid(toprow, gg6, labels = c('', 'C'), label_size = 12, ncol = 1)
dev.off()

################################################## ECOLOGY ANALYSIS ##############################

##### RUN LIFE HISTORY PGLS #####
microdata <- eardata %>% filter(!is.na(Microhabitat))
table(microdata$ear)
microdata <- microdata %>% filter(Microhabitat != 'aquatic-semi.arboreal' &
                                  Microhabitat != 'burrowing-semi.arboreal' &
                                  Microhabitat != 'terrestrial-torrential' & 
                                  Microhabitat != 'aquatic-burrowing')
microdata <- microdata %>% mutate(ear = if_else(ear== "Y", 1, 0))
microdata<-as.data.frame(microdata)
row.names(microdata)<-microdata$Genus_species
microtree<-trimtree(microdata, AnuranTree)
microdata<-trimdata(microdata, microtree)
table(microdata$ear)

nosem_microdata <- microdata %>% filter(Microhabitat != 'semi.aquatic' &
                                    Microhabitat != 'semi.arboreal' & 
                                    Microhabitat != 'semi.burrowing')
nosem_microtree<-trimtree(nosem_microdata, AnuranTree)

activedata <- eardata %>% filter(!is.na(ActiveTime), ActiveTime != '__', ActiveTime != 'D_C_N', ActiveTime != 'D__N', ActiveTime != '_C_')
activedata <- activedata %>% mutate(ear = if_else(ear== "Y", 1, 0))
table(activedata$ear)
row.names(activedata)<-activedata$Genus_species
activedata$ActiveTime <- gsub("D_C_", "Diurnal or Diurnal/Crepuscular", activedata$ActiveTime)
activedata$ActiveTime <- gsub("D__", "Diurnal or Diurnal/Crepuscular", activedata$ActiveTime)
activedata$ActiveTime <- gsub("_C_N", "Crepuscular/Nocturnal or Nocturnal", activedata$ActiveTime)
activedata$ActiveTime <- gsub("__N", "Crepuscular/Nocturnal or Nocturnal", activedata$ActiveTime)
#activedata$ActiveTime <- gsub("_C_", "Crepuscular only", activedata$ActiveTime)
activetree<-trimtree(activedata, AnuranTree)
activedata<-trimdata(activedata, activetree)
table(activedata$ear)
table(activedata$ActiveTime)
table(activedata$ActiveTime)[1]/length(activedata[,1])
table(activedata$ActiveTime)[2]/length(activedata[,1])
table(activedata$ActiveTime)[3]/length(activedata[,1])


microcaperdata<-caperdata(microdata, microtree)
nosem_microcaperdata<-caperdata(nosem_microdata, nosem_microtree)
activecaperdata<-caperdata(activedata, activetree)

micro.pgls<-pgls(ear~Microhabitat, data= microcaperdata, lambda="ML", bounds=list(lambda=c(0.10,0.99)))
nosem_micro.pgls<-pgls(ear~Microhabitat, data= nosem_microcaperdata, lambda="ML", bounds=list(lambda=c(0.10,0.99)))
active.pgls<-pgls(ear~ActiveTime, data= activecaperdata, lambda="ML", bounds=list(lambda=c(0.1,0.99)))

plot(micro.pgls)
plot(nosem_micro.pgls)
plot(active.pgls)

summary(micro.pgls)
summary(nosem_micro.pgls)
summary(active.pgls)

##### PLOT PROPORTION PLOTS #####
#Remove excess data from species not found in tree
microdata <- eardata %>% filter(!is.na(Microhabitat))
table(microdata$ear)
microdata <- microdata %>% filter(Microhabitat != 'aquatic-semi.arboreal' &
                                    Microhabitat != 'burrowing-semi.arboreal' &
                                    Microhabitat != 'terrestrial-torrential' & 
                                    Microhabitat != 'aquatic-burrowing') #remove ambiguous microhabitat states
microdata<-as.data.frame(microdata)
row.names(microdata)<-microdata$Genus_species
#microtree<-trimtree(microdata, AnuranTree)
#microdata<-trimdata(microdata, microtree)
table(microdata$ear)

microprop <- dplyr::select(microdata, c(Genus_species, ear, Microhabitat))
table(microprop$Microhabitat)
microdata<-getprops(microprop, Microhabitat)
microdata$Microhabitat <- factor(microdata$Microhabitat, levels = c("arboreal","semi.arboreal","terrestrial","semi.burrowing","burrowing","torrential","semi.aquatic","aquatic"))
microdata$ear <- factor(microdata$ear, levels = c("Y", "N"))
micropal <- c("forestgreen", "chartreuse2","goldenrod2","pink","red4","purple3","lightblue","dodgerblue2","grey70")

activedata <- eardata %>% filter(!is.na(ActiveTime), ActiveTime != '__')
table(activedata$ear)
row.names(activedata)<-activedata$Genus_species
#activetree<-trimtree(activedata, AnuranTree)
#activedata<-trimdata(activedata, activetree)
table(activedata$ear)
activprop <- dplyr::select(activedata, c(Genus_species, ear, ActiveTime))
table(activprop$ActiveTime)

activedata<-getprops(activprop, ActiveTime)
activedata$ActiveTime <- gsub("D_C_N", "Diurnal/Crepuscular/Nocturnal or Diurnal/Nocturnal", activedata$ActiveTime)
activedata$ActiveTime <- gsub("D__N", "Diurnal/Crepuscular/Nocturnal or Diurnal/Nocturnal", activedata$ActiveTime)
activedata$ActiveTime <- gsub("D_C_", "Diurnal or Diurnal/Crepuscular", activedata$ActiveTime)
activedata$ActiveTime <- gsub("D__", "Diurnal or Diurnal/Crepuscular", activedata$ActiveTime)
activedata$ActiveTime <- gsub("_C_N", "Crepuscular/Nocturnal or Nocturnal", activedata$ActiveTime)
activedata$ActiveTime <- gsub("__N", "Crepuscular/Nocturnal or Nocturnal", activedata$ActiveTime)
activedata$ActiveTime <- gsub("_C_", "Crepuscular", activedata$ActiveTime)
activedata$ActiveTime <- factor(activedata$ActiveTime, levels = c("Diurnal or Diurnal/Crepuscular","Crepuscular","Crepuscular/Nocturnal or Nocturnal", "Diurnal/Crepuscular/Nocturnal or Diurnal/Nocturnal"))
activedata$ear <- factor(activedata$ear, levels = c("Y", "N"))
activepal <- c( "yellow","orangered","black","grey90")
#activedata$ActiveTime <- factor(activedata$ActiveTime, levels = c("DCN or DN","D","DC","C","CN","N"))
#activepal <- c( "grey90","lightblue","skyblue","blue","navy","black", "grey70")

devodata <- eardata %>% filter(!is.na(Concise_DevoMode))
devodata <- devodata %>% filter(Concise_DevoMode != '')
table(devodata$ear)
row.names(devodata)<-devodata$Genus_species
#devotree<-trimtree(devodata, AnuranTree)
#devodata<-trimdata(devodata, devotree)
table(devodata$ear)
devoprop <- dplyr::select(devodata, c(Genus_species, ear, Concise_DevoMode))
table(devoprop$Concise_DevoMode)
devodata<-getprops(devoprop, Concise_DevoMode)
devodata$Concise_DevoMode <- factor(devodata$Concise_DevoMode, levels = c("Viv","Dir","Lar"))
devopal <- c("maroon","red","blue")

propplot <- function(data, x, y, z, colrpal) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{z}})) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colrpal)+
    theme_bw() + 
    theme(legend.position="top", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}

pp1<-propplot(microdata, ear, prop, Microhabitat, micropal)
pp2<-propplot(activedata, ear, prop, ActiveTime, activepal)
pp3<-propplot(devodata, ear, prop, Concise_DevoMode, devopal)

plot_grid(pp1, pp2, pp3, labels = c('A', 'B','C'), ncol=1)

pdf(file = "Fig2_PropPlots_.pdf", width = 12,  height = 6)
plot_grid(pp1, pp2, pp3, labels = c('A','B','C'), ncol=3)
dev.off()



################################################## CALL & SENSORY CONSTRAINTS ####################

##### RUN BODY SIZE AND CALL PGLS #####
#Run call dominant frequency PGLS
calldata <- eardata %>% filter(!is.na(df))
calldata <- eardata %>% filter(df > 0)
calldata <- calldata %>% mutate(ear = if_else(ear== "Y", 1, 0))
#calldata$df <- as.numeric(calldata$df)
length(calldata[,1])
table(calldata$ear)

row.names(calldata)<-calldata$Genus_species
calltree<-trimtree(calldata, AnuranTree)
calldata<-trimdata(calldata, calltree)
table(calldata$ear)

callcaperdata<-caperdata(calldata, calltree)
call.pgls<-pgls(ear~log(df), data= callcaperdata, lambda="ML", bounds=list(lambda=c(0.10,0.99)))

plot(call.pgls)

summary(call.pgls)

#Run body size PGLS
sizedata <- eardata %>% filter(!is.na(Collapsed_SVL))
sizedata <- sizedata %>% filter(Collapsed_SVL>0)
sizedata <- sizedata %>% mutate(ear = if_else(ear== "Y", 1, 0))
length(sizedata[,1])
table(sizedata$ear)

row.names(sizedata)<-sizedata$Genus_species
sizetree<-trimtree(sizedata, AnuranTree)
sizedata<-trimdata(sizedata, sizetree)
table(sizedata$ear)

sizecaperdata<-caperdata(sizedata, sizetree)
size.pgls<-pgls(ear~log(Collapsed_SVL), data= sizecaperdata, lambda="ML", bounds=list(lambda=c(0.30,0.99)))

plot(size.pgls)

summary(size.pgls)

#Run call ~ body callBS PGLS
callBSdata <- sizedata %>% filter(!is.na(df))
length(callBSdata[,1])
table(callBSdata$ear)

callBSdata<-as.data.frame(callBSdata)
row.names(callBSdata)<-callBSdata$Genus_species
callBStree<-trimtree(callBSdata, AnuranTree)
#callBSdata<-trimdata(callBSdata, callBStree)
table(callBSdata$ear)

callBScaperdata<-caperdata(callBSdata, callBStree)
callBS.pgls<-pgls(ear~log(Collapsed_SVL)*log(df), data= callBScaperdata, lambda="ML", bounds=list(lambda=c(0.10,0.99)))
plot(callBS.pgls)
summary(callBS.pgls)
anova(callBS.pgls)

##### PLOT EAR LOSS vs BODY SIZE, EAR LOSS vs CALL DOM FREQ, BODY SIZE vs CALL DOM FREQ w/ EAR LOSS COLOR #####
calldata <- eardata %>% filter(!is.na(df))
calldata <- eardata %>% filter(df > 0)
calldata <- calldata %>% mutate(ear = if_else(ear== "Y", 1, 0))
#calldata$df <- as.numeric(calldata$df)
length(calldata[,1])
table(calldata$ear)

sizedata <- eardata %>% filter(!is.na(Collapsed_SVL))
sizedata <- sizedata %>% filter(Collapsed_SVL>0)
sizedata <- sizedata %>% mutate(ear = if_else(ear== "Y", 1, 0))
length(sizedata[,1])
table(sizedata$ear)

callBSdata <- sizedata %>% filter(!is.na(df))
length(callBSdata[,1])
table(callBSdata$ear)


callplot1 <- function(data, x, y) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{x}})) + 
  geom_violin(alpha = 0.5) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) + 
  theme_classic() + 
  scale_fill_manual(values = pal) 
}

callplot2 <- function(data, x, y, sumstat) {
  ggplot(data, aes(x = {{x}}, y = {{y}}, fill = {{x}})) + 
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .25, outlier.shape = NA) +
    geom_point(size = 1, alpha = .2, position = position_jitter(seed = 1, width = .1)) + 
    annotation_custom(grobTree(textGrob(paste(" F = ", round(sumstat$fstatistic[1], digits = 2), "\n", "df = ", sumstat$fstatistic[2],"/", sumstat$fstatistic[3], "\n", "p = ", round(sumstat$coefficients[2,4], digits = 3)), x=0.4,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=12, fontface="italic")))) +
    theme_classic() + 
    theme(legend.position="none") +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) 
}

sizedata$ear <- as.factor(sizedata$ear)
calldata$ear <- as.factor(calldata$ear)
callBSdata$ear <- as.factor(callBSdata$ear)

callBSdata$ear <- factor(callBSdata$ear, levels = c('1', '0'), ordered = T)
calldata$ear <- factor(calldata$ear, levels = c('1', '0'), ordered = T)
sizedata$ear <- factor(sizedata$ear, levels = c('1', '0'), ordered = T)
pal <- c("#FF8C00", "#2D7297")

gg1<-callplot2(sizedata, ear, log(Collapsed_SVL), summary(size.pgls))
gg2<-callplot2(calldata, ear, log(df), summary(call.pgls))
sumstat<-anova(callBS.pgls)
gg3<-ggplot(callBSdata, aes(x = log(Collapsed_SVL), y = log(df), color = ear)) +
  geom_point(size = 3, alpha = .7, position = position_jitter(seed = 1, width = 0.2)) +
  geom_smooth(method = "lm", fill = NA) +
  annotation_custom(grobTree(textGrob(paste(" F = ", round(sumstat$`F value`[3], digits = 2), "\n", "df = ", sumstat$Df[3],"/", sumstat$Df[4], "\n", "p = ", round(sumstat$`Pr(>F)`[3], digits = 3)), x=0.7,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=12, fontface="italic")))) +
  scale_color_manual(values=pal) +
  theme_classic()

toprow <- plot_grid(gg1, gg2, labels = c('A', 'B'), label_size = 12)
plot_grid(toprow, gg3, labels = c('', 'C'), label_size = 12, ncol = 1)


pdf(file = "Fig3_CallPlots_.pdf", width = 8,  height = 8)
plot_grid(toprow, gg3, labels = c('', 'C'), label_size = 12, ncol = 1)
dev.off()






















