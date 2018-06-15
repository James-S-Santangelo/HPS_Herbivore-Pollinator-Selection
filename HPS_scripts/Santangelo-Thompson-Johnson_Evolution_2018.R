###################################################################
#### HERBIVORES AND PLANT DEFENSES AFFECT SELECTION ON PLANT ######
#### REPRODUCTIVE TRAITS MORE STRONGLY THAN POLLINATORS ###########
#################################################################### 

# Authors: James S. Santangelo, Ken A. Thompson and Marc T. J. Johnson
# Journal: Journal of Evolutionary Biology
# Year: 2018
# Volume:
# Issue:
# Pages:

###############
#### SETUP ####
###############

# Restart R
.rs.restartR()
# Clear workspace
rm(list = ls())

# Create checkpoint with package versions on date analysis was performed.
# Install packages and dependencies in project root.
# Run using R v.3.4.3.
# Checkpoint will ask to create directory. Answer 'y'. 
# Do not proceed with running remaining code until this directory is created.
library(checkpoint)
checkpoint("2018-02-10", project = getwd(),
           checkpointLocation = "./", verbose = TRUE, 
           forceProject = TRUE)

# Confirm that checkpoint worked
getOption("repos") # SHould return MRAN mirror with date = 2018-02-10
normalizePath(.libPaths(), winslash = "/") # Should return HPS project .checkpoint path before default system R library path
installed.packages(.libPaths()[1])[, "Package"] # Show installed packages

#Load required packages
library(Rmisc)
library(lme4)
library(lmerTest)
library(lsmeans)
library(ggplot2)
library(broom)
library(car)
library(multcompView)
library(dplyr)

#FUNCTIONS

#Function for calculating broad-sense heritatility.
Broad_sense <- function(x) {
  x$vcov[1]/sum(x$vcov)
}

##Function for calculating coefficient of genotypic variance
#x is a variance dataframe and y is trait mean
CVg <- function(x,y) {
  100 * (sqrt(x$vcov[1])/y)
}


#Load data on all plants. Re-name columns, remove dead plants, set contrasts, set factors.
Experimental.data = "HPS_data-clean/Santangelo-Thompson-Johnson_Evolution_2018_Experimental-data.txt"
datExp <- read.table(Experimental.data, header = T,na.strings=c("NA", "#DIV/0!"), fill = T)
datExp <- within(datExp, {
  HCN = ifelse(HCN == 0, "No", "Yes")
})
datExp$HCN <- as.factor(datExp$HCN)
datExp$Glycosides.Ac <- as.factor(datExp$Glycosides.Ac)
datExp$Linamarase.Li <- as.factor(datExp$Linamarase.Li)
datExp <- datExp[-which(datExp$Status == "DEAD"),]
names(datExp)[names(datExp) == "Parent.Genotype"] <- "Genotype"
options(contrasts = c("contr.sum","contr.poly"))
datExp <- within(datExp, {
  Mammal.herb = ifelse(is.na(Mammal.herb), 0, 1)
})
datExp$Mammal.herb <- as.factor(datExp$Mammal.herb)
datExp <- within(datExp, {
  Glycosides.Ac = ifelse(Glycosides.Ac == 1, "Yes", "No")
  Linamarase.Li = ifelse(Linamarase.Li == 1, "Yes", "No")
})
datExp$Glycosides.Ac <- as.factor(datExp$Glycosides.Ac)
datExp$Linamarase.Li <- as.factor(datExp$Linamarase.Li)

#Theme used for plots throughout script
ng1=theme(aspect.ratio=0.7,panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line.x = element_line(color="black",size=1),
          axis.line.y = element_line(color="black",size=1),
          axis.ticks=element_line(color="black"),
          axis.text=element_text(color="black",size=15),
          axis.title=element_text(color="black",size=1),
          axis.title.y=element_text(vjust=2,size=17),
          axis.title.x=element_text(vjust=0.1,size=17),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          strip.text.x = element_text(size = 10, colour = "black",face = "bold"),
          strip.background = element_rect(colour="black"),
          legend.position = "top", legend.direction="vertical",
          legend.text=element_text(size=17), legend.key = element_rect(fill = "white"),
          legend.title = element_text(size=17),legend.key.size = unit(1.0, "cm"))

ng1.45=theme(aspect.ratio=0.7,panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border=element_blank(),
             axis.line.x = element_line(color="black",size=1),
             axis.line.y = element_line(color="black",size=1),
             axis.ticks=element_line(color="black"),
             axis.text=element_text(color="black",size=15),
             axis.title=element_text(color="black",size=1),
             axis.title.y=element_text(vjust=2,face="bold",size=18),
             axis.title.x=element_text(vjust=0.1,face="bold",size=18),
             axis.text.x=element_text(size=16,angle=45,hjust=1),
             axis.text.y=element_text(size=16),
             legend.position = "top", legend.direction="vertical",
             legend.text=element_text(size=13), legend.key = element_rect(fill = "white"),
             legend.title = element_text(size=15,face="bold"),legend.key.size = unit(1.0, "cm"))

###############################################################################
#### SUPPLEMENTARY TEXT: EFFECTS OF PESTICIDES ON PLANT GROWTH AND FITNESS ####
###############################################################################

#Load insecticide trial dataset
Insecticide.data = "HPS_data-clean/Santangelo-Thompson-Johnson_Evolution_2018_Insecticide-trial.txt"
datIns <- read.table(Insecticide.data, header=T, fill = T, row.names = NULL)

#Mixed model testing the effects of insecticide and molluscicide on plant biomass
Model.Ins.Biomass <- lmer(Biomass~Molluscicide*Insecticide + (1|Parent.plant), data = datIns)
summary(Model.Ins.Biomass)
anova(Model.Ins.Biomass, ddf = "kenward-roger", type=3)

#Post-hoc comparisons of mixed model using Kenward-Roger degrees of freedom
lsmeans::lsmeans(Model.Ins.Biomass, pairwise~Molluscicide * Insecticide, mode = "kenward-roger")

#Mixed model testing the effects of insecticide and molluscicide on plant fitness
Model.Ins.Fitness <- lmer(Num.Seeds~Molluscicide*Insecticide + (1|Parent.plant), data = datIns, na.action = "na.omit")
summary(Model.Ins.Fitness)
anova(Model.Ins.Fitness, ddf = "kenward-roger", type=3)

#Post-hoc comparisons of mixed model using Kenward-Roger degrees of freedom
lsmeans::lsmeans(Model.Ins.Fitness, pairwise~Molluscicide * Insecticide, mode = "kenward-roger")

#Figure S2A: Biomass by molluscicide and insecticide
meansBIOxMolxIns <- summarySE(datIns, measurevar = "Biomass", groupvars = c("Molluscicide","Insecticide"), na.rm = T)
meansBIOxMolxIns$Insecticide <- factor(meansBIOxMolxIns$Insecticide, levels=c("None", "Low", "High"))
plotBIOxMolxIns<-ggplot(meansBIOxMolxIns,aes(x=Insecticide, y=Biomass,shape=Molluscicide,fill=Molluscicide))+
  geom_errorbar(aes(ymin=Biomass-se,ymax=Biomass+se),position=position_dodge(width=0.5),width=0.3)+
  geom_point(size = 3, position=position_dodge(width = 0.5))+
  xlab("Insecticide")+ylab("Mean biomass (g)")+
  coord_cartesian(ylim=c(2.5, 4)) + scale_y_continuous(breaks=seq(from=2.5,to=4,by=0.25))+
  scale_shape_manual(labels = c("No molluscicide", "Molluscicide"),values=c(21, 24))+
  scale_fill_manual(labels = c("No molluscicide", "Molluscicide"),values=c("white", "black")) +
  ng1+theme(legend.title=element_blank())

#Figure S2B: Number of seeds by insecticide and molluscicide
meansFitxMolxIns <- summarySE(datIns, measurevar = "Num.Seeds", groupvars = c("Molluscicide","Insecticide"), na.rm = T)
meansFitxMolxIns$Insecticide <- factor(meansFitxMolxIns$Insecticide, levels=c("None", "Low", "High"))
PlotFitxMolxIns<-ggplot(meansFitxMolxIns,aes(x=Insecticide, y=Num.Seeds,shape=Molluscicide,fill=Molluscicide))+
  geom_errorbar(aes(ymin=Num.Seeds-se,ymax=Num.Seeds+se),position=position_dodge(width=0.5),width=0.3)+
  geom_point(size = 3, position=position_dodge(width = 0.5))+
  xlab("Insecticide")+ylab("Mean number of seeds")+
  coord_cartesian(ylim=c(20,90))+scale_y_continuous(breaks=seq(from=20,to=90,by=10))+
  scale_shape_manual(labels = c("No molluscicide", "Molluscicide"),values=c(21, 24))+
  scale_fill_manual(labels = c("No molluscicide", "Molluscicide"),values=c("white", "black")) +
  ng1+theme(legend.title=element_blank())

# Save figures 2A and 2B to current working directory
ggsave("HPS_figures/Figure.S2A_Biomass.x.Ins.pdf", plot = plotBIOxMolxIns, width = 5, height = 5, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.S2B_Fitness.x.Ins.pdf", plot = PlotFitxMolxIns, width = 5, height = 5, unit = "in", dpi = 600)


############################################################################
#### SUPPLEMENTARY TEXT: EFFECTS OF PESTICIDES ON POLLINATOR VISITATION ####
############################################################################

#Load pollinator obervation dataset
Pollinator.observations = "HPS_data-clean/Santangelo-Thompson-Johnson_Evolution_2018_Pollinator-observations.txt"
datPoll.obs <- read.table(Pollinator.observations, header = T, fill = T, row.names = NULL)

#LMM testing for effects of insecticide application on time spent foraging on plants. Control for plant display
model.Poll.obs <- lmer(Time.s ~ Num.Inf + Insecticide + (1|Pollinator/Plant) + (1|Date), data = datPoll.obs)
summary(model.Poll.obs)
anova(model.Poll.obs, type = 3, ddf = "kenward-roger")

#Effects of insecticides on visitation
chisq.test(table(datPoll.obs$Insecticide))
datPoll.obs %>% group_by(Insecticide) %>% dplyr::summarize(count = n())

###########################################################
#### SUPPLEMENTART TEXT: EFFICACY OF HAND POLLINATIONS ####
###########################################################

#Create dataset containing only plants that had an inflorescence bagged
datBag <- subset(datExp,Bag.Seed.Num > 0)
datBag$Seed.mass.Inf <- datBag$Seed.only / datBag$Total.Inf

#Melt dataset for analysis
datBag.Melt <- dplyr::select(datBag,Plant,Stolon,Herbivory,HCN,Pollination,Bag.Seed.Wht, Seed.mass.Inf)
datBag.Melt <- data.table::melt(datBag.Melt,id.vars = c("Plant","Pollination","Herbivory",
                                                        "HCN","Stolon"),
                            variable.name = "Bagged",
                            value.name = "Seed.mass")
datBag.Melt <- within(datBag.Melt, {
  Bagged = ifelse(Bagged == "Bag.Seed.Wht","Yes","No")})
datBag.Melt <- within(datBag.Melt, {
  Trt = ifelse(Pollination == "Open" & Bagged == "No", "Open pollinated",
         ifelse(Pollination == "Supp" & Bagged == "No", "Open + hand pollinated", "Hand pollinated only"))
})

#Model seeing if the mass of seeds from bagged inflorescences differs from mass of
#seeds from unbagged inflorescences from open pollinated plants
Mod.Bag <- lmer(Seed.mass ~ Trt + (1|Plant/Bagged), data = datBag.Melt)
anova(Mod.Bag, type = 3, ddf = "kenward-roger")
summary(Mod.Bag)

#Post-hoc comparisons of above model
lsmeans::lsmeans(Mod.Bag, pairwise~Trt)

#Means
meansSeed.x.Trt <- summarySE(datBag.Melt, measurevar = "Seed.mass", groupvars = 'Trt')

#Figure S3: Plot seed mass of bagged and ubagged inflorescences for open and hand pollinated plants
#on the same figure
levels = c("Open pollinated","Hand pollinated only","Open + hand pollinated")
meansSeed.x.Trt$Trt <- factor(meansSeed.x.Trt$Trt, levels = levels)
plotSeed.x.Trt <-ggplot(datBag.Melt, aes(x=Trt, y=Seed.mass))+
  geom_boxplot()+
  # geom_errorbar(aes(ymin=Seed.mass-se,ymax=Seed.mass+se),width=0.15)+
  xlab("")+ylab("Mean seed
mass (g)") + ng1.45
plotSeed.x.Trt

#Save figure S3
ggsave("HPS_figures/Figure.S3_Fitness.x.Poll.x.Bagged.pdf", plot = plotSeed.x.Trt, width = 5, height = 5, unit = "in", dpi = 600)

###########################
#### HERBIVORY SURVEYS ####
###########################

## MIXED MODELS TESTING FOR DIFFERENCES IN DAMAGE ##

#Sample sizes for herbivory survey
datExp %>% group_by(HCN) %>% summarize(n = sum(!is.na(Leaf.avg.dmg.1)))
datExp %>% group_by(Herbivory) %>% summarize(n = sum(!is.na(Leaf.avg.dmg.1)))

#Survey 1
modelHerb.1 <- lmer(Leaf.avg.dmg.1~HCN*Herbivory + (1|Genotype) + (1|Block) + (1|Genotype:Herbivory),data = datExp)
anova(modelHerb.1, type = 3, ddf = 'kenward-roger')
summary(modelHerb.1)
Tukey_modelHerb.1 <- lsmeans::lsmeans(modelHerb.1, pairwise~HCN*Herbivory, mode = "kenward-roger")
lsmeans::cld(Tukey_modelHerb.1)

#Survey 2
modelHerb.2 <- lmer(Leaf.avg.dmg.2~HCN*Herbivory + (1|Genotype) + (1|Block) + (1|Genotype:Herbivory),data = datExp)
anova(modelHerb.2, type = 3, ddf = 'kenward-roger')
summary(modelHerb.2)
Tukey_modelHerb.2 <- lsmeans::lsmeans(modelHerb.2, pairwise~HCN*Herbivory, mode = "kenward-roger")
lsmeans::cld(Tukey_modelHerb.2)

#Survey 3
modelHerb.3 <- lmer(Leaf.avg.dmg.3~HCN*Herbivory + (1|Genotype) + (1|Block) + (1|Genotype:Herbivory),data = datExp)
anova(modelHerb.3, type = 3, ddf = 'kenward-roger')
summary(modelHerb.3)
Tukey_modelHerb.3 <- lsmeans::lsmeans(modelHerb.3, pairwise~HCN*Herbivory, mode = "kenward-roger")
lsmeans::cld(Tukey_modelHerb.3)

#3 surveys combined
datExp$Leaf.avg.dmg.All <- (datExp$Leaf.avg.dmg.1 + datExp$Leaf.avg.dmg.2 + datExp$Leaf.avg.dmg.3)/3
modelHerb.4 <- lmer(Leaf.avg.dmg.All~HCN*Herbivory + (1|Genotype) + (1|Genotype:Herbivory),data = datExp)
anova(modelHerb.4, type = 3, ddf = 'kenward-roger')
summary(modelHerb.4)

#Create melted dataset for plotting all 3 surveys on single figure
datHerb <- dplyr::select(datExp,Genotype,Stolon,Herbivory,HCN,Leaf.avg.dmg.1,
  Leaf.avg.dmg.2,Leaf.avg.dmg.3)
datHerb <- data.table::melt(datHerb,id.vars = c("Genotype","Stolon","Herbivory","HCN"),
  variable.name = "Survey",
  value.name = "Damage")
datHerb <- within(datHerb, {
  Survey = ifelse(Survey == "Leaf.avg.dmg.1","Early",
            ifelse(Survey == "Leaf.avg.dmg.2","Mid","Late"))
})

#Mean herbivore damage across all surveys with SE
mean(datHerb$Damage, na.rm = T)
mean(datHerb$Damage, na.rm = T)/sqrt(sd(datHerb$Damage, na.rm = T))


#Means for herbivory surveys
Summary.Herb.x.Herb_HCN <- summarySE(datHerb, measurevar = "Damage", groupvar = c("HCN", "Herbivory","Survey"), na.rm = T)
Summary.Herb.x.HCN <- summarySE(datHerb, measurevar = "Damage", groupvar = c("HCN","Survey"), na.rm = T)

#Effect sizes (% change). Values from means datasets above
(14.745105 - 5.635190)/14.745105 #Ins, Early
(15.110606 - 4.651139)/15.110606 # Ins, Mid
(9.159596 - 3.667677)/9.159596 # Ins, Late
(10.883080 - 9.500715)/10.883080 # HCN, Early
(11.370633 - 8.408081)/11.370633 # HCN, Late
(7.563038 - 5.270025)/7.563038 # HCN, Mid

#Figure S4. Plot mean damage across HCN and Herbivory for each survey. 
Summary.Herb.x.Herb_HCN$Survey <- as.character(Summary.Herb.x.Herb_HCN$Survey)
Summary.Herb.x.Herb_HCN$Survey <- factor(Summary.Herb.x.Herb_HCN$Survey, levels=c("Early","Mid","Late"))
plot.Herb.x.Herb_HCN <- ggplot(Summary.Herb.x.Herb_HCN,aes(x = HCN, y = Damage,shape = Herbivory,fill = Herbivory, group = Herbivory))+
  geom_errorbar(aes(ymin=Damage-se,ymax=Damage+se),width=0.15,size=0.7)+
  geom_line(size = 1, aes(linetype = Herbivory)) + 
  geom_point(size = 4.5)+
  facet_wrap(  ~ Survey) + 
  xlab("HCN")+ylab("% Herbivore damage")+
  coord_cartesian(ylim = c(1.5,20))+
  scale_y_continuous(breaks = seq(from = 2, to = 20, by = 2))+
  scale_shape_manual(labels = c("Ambient herbivory", "Reduced herbivory"),values=c(22, 23))+
  scale_fill_manual(labels = c("Ambient herbivory", "Reduced herbivory"),values=c("white", "black")) +
  ng1 + theme(aspect.ratio=1.0, legend.title=element_blank())
plot.Herb.x.Herb_HCN

# Save figures S4 to current working directory
ggsave("HPS_figures/Figure.S4_Herbivory.x.Insecticide-HCN_Three.Surveys.pdf", plot = plot.Herb.x.Herb_HCN, width = 10, height = 8, unit = "in", dpi = 600)

#######################
#### FLORAL DAMAGE ####
#######################

#Sample sizes for banner damage
datExp %>% group_by(HCN) %>% summarize(n = sum(!is.na(Avg.bnr.dmg)))
datExp %>% group_by(Herbivory) %>% summarize(n = sum(!is.na(Avg.bnr.dmg)))

#Model testing for differences in floral damage due to HCN and pesticide treatment
model.bnr.dmg <- lmer(Avg.bnr.dmg ~ HCN*Herbivory + (1|Genotype) + (1|Block) + (1|Genotype:Herbivory),data = datExp)
anova(model.bnr.dmg, type = 3, ddf = "kenward-roger")

#Mean floral damage for pesticide treated and untreated plants
means.bnr.dmg <- summarySE(datExp, measurevar = "Avg.bnr.dmg", groupvars = "Herbivory", na.rm = TRUE)
(41.58122 - 14.90667)/41.58122 #ES

#Figure S5A: Banner petal damage for pesticide treated and untreated plants
means.bnr.dmg$Herbivory <- factor(means.bnr.dmg$Herbivory, levels=c("Reduced", "Ambient"))
plot.Bnr.dmg.x.Herb <- ggplot(means.bnr.dmg,aes(x = Herbivory, y = Avg.bnr.dmg))+
  geom_errorbar(aes(ymin=Avg.bnr.dmg-se,ymax=Avg.bnr.dmg+se),width=0.15,size=0.7)+
  geom_point(size = 4.5)+
  xlab("Herbivory")+ylab("% Banner petal damage") +
  coord_cartesian(ylim = c(12, 43.5)) + scale_y_continuous(breaks = seq(from = 15, to = 40, by = 5)) +
  ng1 + theme(legend.title=element_blank())
plot.Bnr.dmg.x.Herb

#Save figure S5
ggsave("HPS_figures/Figure.S5_Bnr.dmg.x.Herbivory.pdf", plot = plot.Bnr.dmg.x.Herb, width = 5, height = 5, unit = "in", dpi = 600)

#####################
#### VOLE DAMAGE ####
#####################

#Total proportion on plants damaged by voles.
table(datExp$Mammal.herb)[2] / sum(table(datExp$Mammal.herb)[1], table(datExp$Mammal.herb)[2])


########################################################
#### TRAIT VARIANCES AND BROAD-SENSE HERITABILITIES ####
########################################################

# Calculation of heritabilities and coefficients of genetic variation for traits in
# our experiment. Done in ambient and suppressed herbivore treatments separatly. Not
# done in pollination treatments separately since genotypes did not vary significantly
# in response to pollination treatment (see trait models in next section).
# Values are those in Table S3

datExp_Amb <- subset(datExp, Herbivory == "Ambient")
datExp_Red <- subset(datExp, Herbivory == "Reduced")

## AMBIENT HERBIVORY

##Calculate significant of genotype (i.e. Plant)
#Run mixed modes with genotype and block. REML since unbiased for random effects
model.1.full_Amb <- lmer(Flower.date ~ (1|Block) + (1|Genotype), data = datExp_Amb, REML = T)
model.2.full_Amb <- lmer(Avg.Bnr.Wdth ~ (1|Block) + (1|Genotype), data = datExp_Amb, REML = T)
model.3.full_Amb <- lmer(Avg.Bnr.Ht ~ (1|Block) + (1|Genotype), data = datExp_Amb, REML = T)
model.4.full_Amb <- lmer(Total.Inf ~ (1|Block) + (1|Genotype), data = datExp_Amb, REML = T)
model.5.full_Amb <- lmer(Biomass.plant ~ (1|Block) + (1|Genotype), data = datExp_Amb, REML = T)
model.6.full_Amb <- lmer(Num.flwrs ~ (1|Block) + (1|Genotype), data = datExp_Amb, REML = T)
model.7.full_Amb <- lmer(Total.Seed.mass ~ (1|Block) + (1|Genotype), data = datExp_Amb, REML = T)

#Likelihood ratio tests for all random effect
rand(model.1.full_Amb) # Flower date
rand(model.2.full_Amb) # Banner width
rand(model.3.full_Amb) # Banner height
rand(model.4.full_Amb) # Total inflorescences
rand(model.5.full_Amb) # Vegetative biomass
rand(model.6.full_Amb) # Number of flowers
rand(model.7.full_Amb) # Total seed mass

##Create dataset from random effect components of full models
#Used to estimate genetic variance, heritability and CVg
Var.model.1_Amb <- as.data.frame(VarCorr(model.1.full_Amb),comp="Variance")
Var.model.2_Amb <- as.data.frame(VarCorr(model.2.full_Amb),comp="Variance")
Var.model.3_Amb <- as.data.frame(VarCorr(model.3.full_Amb),comp="Variance")
Var.model.4_Amb <- as.data.frame(VarCorr(model.4.full_Amb),comp="Variance")
Var.model.5_Amb <- as.data.frame(VarCorr(model.5.full_Amb),comp="Variance")
Var.model.6_Amb <- as.data.frame(VarCorr(model.6.full_Amb),comp="Variance")
Var.model.7_Amb <- as.data.frame(VarCorr(model.7.full_Amb),comp="Variance")

#Broad sense heritabilities
Broad_sense(Var.model.1_Amb) # Flowering date
Broad_sense(Var.model.2_Amb) # Banner width
Broad_sense(Var.model.3_Amb) # Banner length
Broad_sense(Var.model.4_Amb) # Number of inflorescences
Broad_sense(Var.model.5_Amb) # Biomass
Broad_sense(Var.model.6_Amb) # Number of flowers
Broad_sense(Var.model.7_Amb) # Seed mass

#Trait means
mean.1_Amb <- mean(datExp_Amb$Flower.date, na.rm = T)
mean.2_Amb <- mean(datExp_Amb$Avg.Bnr.Wdth, na.rm = T)
mean.3_Amb <- mean(datExp_Amb$Avg.Bnr.Ht, na.rm = T)
mean.4_Amb <- mean(datExp_Amb$Total.Inf, na.rm = T)
mean.5_Amb <- mean(datExp_Amb$Biomass.plant, na.rm = T)
mean.6_Amb <- mean(datExp_Amb$Num.flwrs, na.rm = T)
mean.7_Amb <- mean(datExp_Amb$Total.Seed.mass, na.rm = T)

#Coefficients of genotypic variance
CVg(Var.model.1_Amb,mean.1_Amb) # Flowering date
CVg(Var.model.2_Amb,mean.2_Amb) # Banner width
CVg(Var.model.3_Amb,mean.3_Amb) # Banner height
CVg(Var.model.4_Amb,mean.4_Amb) # Number of inflorescences
CVg(Var.model.5_Amb,mean.5_Amb) # Biomass
CVg(Var.model.6_Amb,mean.6_Amb) # Number of flowers
CVg(Var.model.7_Amb,mean.7_Amb) # Seed mass

## REDUCED HERBIVORY

##Calculate significant of genotype (i.e. Plant)
#Run mixed modes with genotype and block. REML since unbiased for random effects
model.1.full_Red <- lmer(Flower.date ~ (1|Block) + (1|Genotype), data = datExp_Red, REML = T)
model.2.full_Red <- lmer(Avg.Bnr.Wdth ~ (1|Block) + (1|Genotype), data = datExp_Red, REML = T)
model.3.full_Red <- lmer(Avg.Bnr.Ht ~ (1|Block) + (1|Genotype), data = datExp_Red, REML = T)
model.4.full_Red <- lmer(Total.Inf ~ (1|Block) + (1|Genotype), data = datExp_Red, REML = T)
model.5.full_Red <- lmer(Biomass.plant ~ (1|Block) + (1|Genotype), data = datExp_Red, REML = T)
model.6.full_Red <- lmer(Num.flwrs ~ (1|Block) + (1|Genotype), data = datExp_Red, REML = T)
model.7.full_Red <- lmer(Total.Seed.mass ~ (1|Block) + (1|Genotype), data = datExp_Red, REML = T)

#Likelihood ratio tests for all random effect
rand(model.1.full_Red) # Flower date
rand(model.2.full_Red) # Banner width
rand(model.3.full_Red) # Banner height
rand(model.4.full_Red) # Total inflorescences
rand(model.5.full_Red) # Vegetative biomass
rand(model.6.full_Red) # Number of flowers
rand(model.7.full_Red) # Total seed mass

##Create dataset from random effect components of full models
#Used to estimate genetic variance, heritability and CVg
Var.model.1_Red <- as.data.frame(VarCorr(model.1.full_Red),comp="Variance")
Var.model.2_Red <- as.data.frame(VarCorr(model.2.full_Red),comp="Variance")
Var.model.3_Red <- as.data.frame(VarCorr(model.3.full_Red),comp="Variance")
Var.model.4_Red <- as.data.frame(VarCorr(model.4.full_Red),comp="Variance")
Var.model.5_Red <- as.data.frame(VarCorr(model.5.full_Red),comp="Variance")
Var.model.6_Red <- as.data.frame(VarCorr(model.6.full_Red),comp="Variance")
Var.model.7_Red <- as.data.frame(VarCorr(model.7.full_Red),comp="Variance")

#Broad sense heritabilities
Broad_sense(Var.model.1_Red) # Flowering date
Broad_sense(Var.model.2_Red) # Banner width
Broad_sense(Var.model.3_Red) # Banner length
Broad_sense(Var.model.4_Red) # Number of inflorescences
Broad_sense(Var.model.5_Red) # Biomass
Broad_sense(Var.model.6_Red) # Number of flowers
Broad_sense(Var.model.7_Red) # Seed mass

#Trait means
mean.1_Red <- mean(datExp_Red$Flower.date, na.rm = T)
mean.2_Red <- mean(datExp_Red$Avg.Bnr.Wdth, na.rm = T)
mean.3_Red <- mean(datExp_Red$Avg.Bnr.Ht, na.rm = T)
mean.4_Red <- mean(datExp_Red$Total.Inf, na.rm = T)
mean.5_Red <- mean(datExp_Red$Biomass.plant, na.rm = T)
mean.6_Red <- mean(datExp_Red$Num.flwrs, na.rm = T)
mean.7_Red <- mean(datExp_Red$Total.Seed.mass, na.rm = T)

#Coefficients of genotypic variance
CVg(Var.model.1_Red,mean.1_Red) # Flowering date
CVg(Var.model.2_Red,mean.2_Red) # Banner width
CVg(Var.model.3_Red,mean.3_Red) # Banner height
CVg(Var.model.4_Red,mean.4_Red) # Number of inflorescences
CVg(Var.model.5_Red,mean.5_Red) # Biomass
CVg(Var.model.6_Red,mean.6_Red) # Number of flowers
CVg(Var.model.7_Red,mean.7_Red) # Seed mass

########################################
#### EFFECT OF TREATMENTS ON TRAITS ####
########################################

#Note that only final models (post-backward model selection) with transformed
#traits (if application) are shown

#Phenotypic Trait transformations (if applicable) for trait models
datExp$Flower.date.T <- (datExp$Flower.date)^(1/4)
datExp$Biomass.plant.T <- (datExp$Biomass.plant)^(1/4)
datExp$Total.Inf.T <- log1p(datExp$Total.Inf)
datExp$Total.Seed.mass.T <- (datExp$Total.Seed.mass)^(1/3)
datExp$Seeds.Inf.T <- (datExp$Seeds.Inf)^(1/3)

#Create data subsets, removing rows with missing data for traits. Not really necessary.
datExp_Flower.date <- datExp[which(is.finite(datExp$Flower.date)),]
datExp_Avg.Bnr.Wdth <- datExp[which(is.finite(datExp$Avg.Bnr.Wdth)),]
datExp_Avg.Bnr.Ht <- datExp[which(is.finite(datExp$Avg.Bnr.Ht)),]
datExp_Biomass.plant <- datExp[which(is.finite(datExp$Biomass.plant)),]
datExp_Total.Inf <- datExp[which(is.finite(datExp$Total.Inf)),]
datExp_Num.flwrs <- datExp[which(is.finite(datExp$Num.flwrs)),]
datExp_Total.Seed.mass <- datExp[which(is.finite(datExp$Total.Seed.mass)),]

## FLOWERING DATE

#Model for effects of treatments on flowering date
trait.model.1.T.Final <- lmer(Flower.date.T ~ Mammal.herb + Herbivory + (1 | Genotype),
  data = datExp_Flower.date, REML = F)

#Model output
anova(trait.model.1.T.Final, type = 3, ddf = "kenward-roger")
rand(trait.model.1.T.Final) # Significant genotypic variation

#Dataset for plotting effect of voles on flowering date
means.Flwr.date.Vole <- summarySE(datExp, measurevar = "Flower.date", groupvars = "Mammal.herb", na.rm = T)
names(means.Flwr.date.Vole)[names(means.Flwr.date.Vole) == "Mammal.herb"] <- "Voles"
means.Flwr.date.Vole <- within(means.Flwr.date.Vole, {
  Voles = ifelse(Voles == 0, "Undamaged", "Damaged")
})
(41.93182 - 28.94444) # ES

#Dataset for plotting effect of invertebrate herbivores on flowering date. No Plot.
means.Flwr.date.Herb <- summarySE(datExp, measurevar = "Flower.date", groupvars = "Herbivory", na.rm = T)
36.71031 - 30.84570 # ES

#Figure 2D. Effects of voles on flowering date
means.Flwr.date.Vole$Voles <- factor(means.Flwr.date.Vole$Voles, levels=c("Undamaged", "Damaged"))
plot.Flwr.date.x.Voles <- ggplot(means.Flwr.date.Vole,aes(x=Voles, y=Flower.date)) +
  geom_errorbar(aes(ymin=Flower.date-se,ymax=Flower.date+se),width=0.15,size=0.7) +
  geom_point(size=5.5) +
  xlab("Vole damage") + ylab("Days to first flower") +
  coord_cartesian(ylim = c(27,45)) + scale_y_continuous(breaks = seq(from = 27, to = 45, by = 4)) + 
  ng1 + theme(axis.text.x=element_text(size=17), axis.text.y=element_text(size=17))

# Save Figure 2D
ggsave("HPS_figures/Figure.2D_FF.x.Voles.pdf", plot = plot.Flwr.date.x.Voles, width = 5, height = 5, unit = "in", dpi = 600)

## BANNER WIDTH

#Model for effects of treatments on banner width
trait.model.2.UnT.Final <- lmer(Avg.Bnr.Wdth ~ HCN + Herbivory + Pollination + (1 | Genotype) + (1 | Block) + HCN:Pollination,
  data = datExp_Avg.Bnr.Wdth, REML = F)

#Model output
anova(trait.model.2.UnT.Final, type = 3, ddf = "kenward-roger")
rand(trait.model.2.UnT.Final) # Signficant genotypic, block and GT:Herb variation

#Post-hoc comparisons of mixed model using Kenward-Roger degrees of freedom
Tukey_BW.HCN.Poll <- lsmeans::lsmeans(trait.model.2.UnT.Final, pairwise~HCN*Pollination, mode = "kenward-roger")
Tukey_BW.HCN.Poll

#Dataset for plotting effect of invertebrate herbivores Banner width
means.Bnr.wdth.Herb <- summarySE(datExp, measurevar = "Avg.Bnr.Wdth", groupvars = "Herbivory", na.rm = T)
(3.251672 - 3.214632)/3.251672 # ES

#Dataset for plotting effect of HCN and pollinators Banner width
means.Bnr.wdth.HCN.Poll <- summarySE(datExp, measurevar = "Avg.Bnr.Wdth", groupvars = c("HCN","Pollination"), na.rm = T)

#Figure S6: Effects of herbivores on Banner width
plot.Bnr.wdth.x.Herb <- ggplot(means.Bnr.wdth.Herb,aes(x=Herbivory, y=Avg.Bnr.Wdth))+
  geom_errorbar(aes(ymin=Avg.Bnr.Wdth-se,ymax=Avg.Bnr.Wdth+se),width=0.15,size=0.7)+
  geom_point(size=4.5)+
  xlab("Herbivory")+ylab("Banner width (mm)") +
  coord_cartesian(ylim = c(3.185,3.275)) + scale_y_continuous(breaks = seq(from = 3.19, to = 3.27, by = 0.02)) +
  ng1

#Figure S7A: Effects of HCN and pollinators on banner width
plot.Bnr.wdth.x.HCN.Poll <- ggplot(means.Bnr.wdth.HCN.Poll,aes(x=HCN, y=Avg.Bnr.Wdth, shape = Pollination, fill = Pollination))+
  geom_errorbar(aes(ymin=Avg.Bnr.Wdth-se,ymax=Avg.Bnr.Wdth+se),width=0.15,size=0.7,position = position_dodge(width = 0.5))+
  geom_point(size=4.5,position = position_dodge(width = 0.5))+
  xlab("HCN")+ylab("Banner width (mm)")+
  scale_shape_manual(labels = c("Open pollination", "Supplemental pollination"),values=c(21, 24))+
  scale_fill_manual(labels = c("Open pollination", "Supplemental pollination"),values=c("white", "black")) +
  scale_x_discrete(breaks=c("No","Yes"),labels=c("HCN-", "HCN+")) +
  coord_cartesian(ylim = c(3.16,3.34)) + scale_y_continuous(breaks = seq(from = 3.16, to = 3.34, by = 0.02)) +
  ng1 + theme(legend.title=element_blank(),
                                       axis.title.x = element_blank(),
                                       axis.text.x = element_text(face = "bold"))

#Save figures S6 and S7A
ggsave("HPS_figures/Figure.S6_BW.x.Herbivory.pdf", plot = plot.Bnr.wdth.x.Herb, width = 5, height = 5, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.S7A_BW.x.Poll.x.HCN.pdf", plot = plot.Bnr.wdth.x.HCN.Poll, width = 5, height = 5, unit = "in", dpi = 600)

## BANNER LENGTH

#Model for effects of treatments on banner length
trait.model.3.UnT.Final <- lmer(Avg.Bnr.Ht ~ HCN + Herbivory + Pollination +
  (1 | Genotype) + (1 | Block) + (1 | Herbivory:Genotype) +
  HCN:Pollination, data = datExp_Avg.Bnr.Ht, REML = F)

#Model output
anova(trait.model.3.UnT.Final, type = 3, ddf = "kenward-roger")
rand(trait.model.3.UnT.Final) # Signficant genotypic, block and GT:Herb variation

#Post-hoc comparisons of mixed model using Kenward-Roger degrees of freedom
Tukey_BL.HCN.Poll <- lsmeans::lsmeans(trait.model.3.UnT.Final, pairwise~HCN*Pollination, mode = "kenward-roger")
Tukey_BL.HCN.Poll

#Dataset for plotting effect of invertebrate herbivores Banner length
means.Bnr.ht.Herb <- summarySE(datExp, measurevar = "Avg.Bnr.Ht", groupvars = "Herbivory", na.rm = T)
(6.265449 - 6.187059)/6.265449 # ES

#Dataset for plotting effect of HCN and pollinators Banner length
means.Bnr.ht.HCN.Poll <- summarySE(datExp, measurevar = "Avg.Bnr.Ht", groupvars = c("Pollination","HCN"), na.rm = T)

#Figure 2A: Effects of herbivores on Banner length
means.Bnr.ht.Herb$Herbivory <- factor(means.Bnr.ht.Herb$Herbivory, levels=c("Reduced", "Ambient"))
plot.Bnr.ht.x.Herb <- ggplot(means.Bnr.ht.Herb,aes(x=Herbivory, y=Avg.Bnr.Ht))+
  geom_errorbar(aes(ymin=Avg.Bnr.Ht-se,ymax=Avg.Bnr.Ht+se),width=0.15,size=0.7)+
  geom_point(size=5.5)+
  xlab("Herbivory")+ylab("Banner length (mm)") +
  ng1 + theme(axis.text.x=element_text(size=17), axis.text.y=element_text(size=17))

#Figure S7B: Effects of HCN and Pollinators on banner length
plot.Bnr.ht.x.HCN.Poll <- ggplot(means.Bnr.ht.HCN.Poll,aes(x=HCN, y=Avg.Bnr.Ht, shape = Pollination, fill = Pollination))+
  geom_errorbar(aes(ymin=Avg.Bnr.Ht-se,ymax=Avg.Bnr.Ht+se),width=0.15,size=0.7,position = position_dodge(width = 0.5))+
  geom_point(size=4.5,position = position_dodge(width = 0.5))+
  xlab("HCN")+ylab("Banner length (mm)")+
  scale_shape_manual(labels = c("Open pollination", "Supplemental pollination"),values=c(21, 24))+
  scale_fill_manual(labels = c("Open pollination", "Supplemental pollination"),values=c("white", "black")) +
  scale_x_discrete(breaks=c("No","Yes"),labels=c("HCN-", "HCN+")) +
  coord_cartesian(ylim = c(6.10,6.32)) + scale_y_continuous(breaks = seq(from = 6.10, to = 6.32, by = 0.02)) +
  ng1 + theme(legend.title=element_blank(),
                                     axis.title.x = element_blank(),
                                     axis.text.x = element_text(face = "bold"))

#Save figures 2A and S7B
ggsave("HPS_figures/Figure.2A_BW.x.Herbivory.pdf", plot = plot.Bnr.ht.x.Herb, width = 5, height = 5, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.S7B_BW.x.Poll.x.HCN.pdf", plot = plot.Bnr.ht.x.HCN.Poll, width = 5, height = 5, unit = "in", dpi = 600)

## BIOMASS

#Model for effects of treatments on biomass
trait.model.4.T.Final <- lmer(Biomass.plant.T ~ Mammal.herb + HCN + Herbivory +
  Pollination + (1 | Genotype) + (1 | Block) + Mammal.herb:Herbivory +
  HCN:Herbivory + HCN:Pollination + Herbivory:Pollination +
  HCN:Herbivory:Pollination, data = datExp_Biomass.plant, REML = F)

#Model output
anova(trait.model.4.T.Final,type = 3, ddf = "kenward-roger")
rand(trait.model.4.T.Final) # Signficant genotypic and block

#Post-hoc comparisons of mixed model using Kenward-Roger degrees of freedom
Tukey_Bio.HCN.Poll.Herb <- lsmeans::lsmeans(trait.model.4.T.Final, pairwise~HCN*Herbivory*Pollination, mode = "kenward-roger")
Tukey_Bio.HCN.Poll.Herb

#Variation explained by significant terms from above model
trait.model.4.T.Var <- lmer(Biomass.plant.T ~ (1|Mammal.herb) + (1|HCN) + (1|Herbivory) +
  (1|Pollination) + (1 | Genotype) + (1 | Block) + (1|Mammal.herb:Herbivory) +
  (1|HCN:Herbivory) + (1|HCN:Pollination) + (1|Herbivory:Pollination) +
  (1|HCN:Herbivory:Pollination), data = datExp_Biomass.plant, REML = F)
summary(trait.model.4.T.Var)
Biomass.Var <- as.data.frame(VarCorr(trait.model.4.T.Var),comp="Variance")
Biomass.Var[which(Biomass.Var$grp == "HCN:Herbivory:Pollination"),"vcov"]/sum(Biomass.Var$vcov)
Biomass.Var[which(Biomass.Var$grp == "Mammal.herb"),"vcov"]/sum(Biomass.Var$vcov)
Biomass.Var[which(Biomass.Var$grp == "Mammal.herb:Herbivory"),"vcov"]/sum(Biomass.Var$vcov)

#Dataset for plotting effect of vole damage on biomass
means.Biomass.Vole <- summarySE(datExp, measurevar = "Biomass.plant", groupvars = c("Mammal.herb"), na.rm = T)
names(means.Biomass.Vole)[names(means.Biomass.Vole) == "Mammal.herb"] <- "Voles"
means.Biomass.Vole <- within(means.Biomass.Vole, {
  Voles = ifelse(Voles == 0, "Undamaged", "Damaged")
})
(44.78854 - 24.92469)/44.78854 # ES

#Dataset for plotting effect of HCN, herbivores, and pollinators on biomass
means.Biomass.HCN.Herb.Poll <- summarySE(datExp, measurevar = "Biomass.plant", groupvars = c("Herbivory","HCN","Pollination"), na.rm = T)

# Figure S8: Effects of HCN, Herbivores, and pollinators on vegetative biomass
means.Biomass.HCN.Herb.Poll$Pollination <- factor(means.Biomass.HCN.Herb.Poll$Pollination, levels=c("Open", "Supp"))
means.Biomass.HCN.Herb.Poll$Herbivory <- factor(means.Biomass.HCN.Herb.Poll$Herbivory, levels=c("Reduced", "Ambient"))
plot.Biomass.x.Poll.Herb.HCN <- ggplot(means.Biomass.HCN.Herb.Poll,aes(x=Herbivory, y=Biomass.plant, shape = HCN, fill = HCN))+
  geom_errorbar(aes(ymin=Biomass.plant-se,ymax=Biomass.plant+se),width=0.15,size=0.7, position = position_dodge(width = 0.5))+
  geom_point(size=4.5, position = position_dodge(width = 0.5))+
  facet_wrap(  ~ Pollination)+
  xlab("Herbivory")+ylab("Biomass (mg)")+
  scale_shape_manual(labels = c("HCN-","HCN+"),values=c(24, 21), guide = guide_legend(title.position = "top"))+
  scale_fill_manual(labels = c("HCN-","HCN+"),values=c("white", "black")) +
  coord_cartesian(ylim = c(25,52)) + scale_y_continuous(breaks = seq(from = 25, to = 50, by = 5)) +
  ng1 + theme(legend.title = element_blank(),
              legend.direction = "horizontal",
              legend.position = "top",
              panel.spacing = unit(2, "lines"))

#Figure S9. Effects of vole damage on vegetative biomass
plot.Biomass.x.Vole <- ggplot(means.Biomass.Vole,aes(x=Voles, y=Biomass.plant))+
  geom_errorbar(aes(ymin=Biomass.plant-se,ymax=Biomass.plant+se),width=0.15,size=0.7)+
  geom_point(size=4.5)+
  xlab("Vole damage")+ylab("Biomass (mg)") +
  coord_cartesian(ylim = c(20,50)) + scale_y_continuous(breaks = seq(from = 20, to = 50, by = 5)) +
  ng1

#Save figures S8 and S9
ggsave("HPS_figures/Figure.S8_Biomass.x.HCN.Poll.Herb.pdf", plot = plot.Biomass.x.Poll.Herb.HCN, width = 5, height = 5, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.S9_Biomass.x.Voles.pdf", plot = plot.Biomass.x.Vole, width = 5, height = 5, unit = "in", dpi = 600)

## NUMBER OF INFLORESCENCES

#Model for effects of treatments on number of inflorescences
trait.model.5.T.Final <- lmer(Total.Inf.T ~ Mammal.herb + HCN + Herbivory +
  (1 | Genotype) + (1 | Block) + (1 | Herbivory:Genotype) +
  Mammal.herb:HCN, data = datExp_Total.Inf, REML = F)

#Model output
anova(trait.model.5.T.Final, type = 3, ddf = "kenward-roger")
rand(trait.model.5.T.Final) # Significant Genotype, block and GT:Herb variation

#Post-hoc comparisons of mixed model using Kenward-Roger degrees of freedom
Tukey_InfVoles <- lsmeans::lsmeans(trait.model.5.T.Final, pairwise~Mammal.herb*HCN, mode = "kenward-roger")
Tukey_InfVoles

#Dataset for plotting the effects of herbivores on number of inflorescences
means.Inflor.Herb <- summarySE(datExp, measurevar = "Total.Inf", groupvars = "Herbivory", na.rm = T)
(47.52778 - 34.32071)/47.52778 # ES

#Dataset for plotting effects of HCN and Voles on number of inflorescence
means.Inflor.HCN.Voles <- summarySE(datExp, measurevar = "Total.Inf", groupvars = c("HCN","Mammal.herb"), na.rm = T)
names(means.Inflor.HCN.Voles)[names(means.Inflor.HCN.Voles) == "Mammal.herb"] <- "Voles"
means.Inflor.HCN.Voles <- within(means.Inflor.HCN.Voles, {
  Voles = ifelse(Voles == 0, "Undamaged", "Damaged")
})

#Dataset fo getting effects size for effect of vole damage on number of inflorescences
means.Inflor.Voles <- summarySE(datExp, measurevar = "Total.Inf", groupvars = "Mammal.herb", na.rm = T)
(52.21656 - 24.35514)/52.21656

#Figure 2B. Effects of herbivores on number of inflorescences
means.Inflor.Herb$Herbivory <- factor(means.Inflor.Herb$Herbivory, levels=c("Reduced", "Ambient"))
plot.Inflor.x.Herb <- ggplot(means.Inflor.Herb,aes(x=Herbivory, y=Total.Inf))+
  geom_errorbar(aes(ymin=Total.Inf-se,ymax=Total.Inf+se),width=0.15,size=0.7)+
  geom_point(size=5.5)+
  xlab("Herbivory")+ylab("Number of Inflorescences") +
  coord_cartesian(ylim = c(30, 55)) + scale_y_continuous(breaks = seq(from=30, to=55, by=5)) +
  ng1 + theme(axis.text.x=element_text(size=17), axis.text.y=element_text(size=17))

#Figure 2E: Effects of HCN and voles on number of inflorescences
means.Inflor.HCN.Voles$Voles <- factor(means.Inflor.HCN.Voles$Voles, levels=c("Undamaged", "Damaged"))
plot.Inflor.x.Vole.HCN <- ggplot(means.Inflor.HCN.Voles,aes(x=Voles, y=Total.Inf, shape = HCN, fill = HCN))+
  geom_errorbar(aes(ymin=Total.Inf-se,ymax=Total.Inf+se),width=0.15,size=0.7,position = position_dodge(width = 0.5))+
  geom_point(size=5.5, position = position_dodge(width = 0.5))+
  xlab("Voles")+ylab("Number of inflorescences")+
  scale_shape_manual(labels = c("HCN-","HCN+"),values=c(24, 21), guide = guide_legend(title.position = "top"))+
  scale_fill_manual(labels = c("HCN-","HCN+"),values=c("white", "black")) +
  ng1 + theme(legend.title = element_blank()) + theme(axis.text.x=element_text(size=17), axis.text.y=element_text(size=17))

#Save figures 2B and 2E
ggsave("HPS_figures/Figure.2B_Inflorescences.x.Herb.pdf", plot = plot.Inflor.x.Herb, width = 5, height = 5, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.2E_Inflorescences.x.Voles.x.HCN.pdf", plot = plot.Inflor.x.Vole.HCN, width = 5, height = 5, unit = "in", dpi = 600)


## NUMBER OF FLOWERS

#Model for effects of treatments on number of flowers per inflorescence
trait.model.6.UnT.Final <- lmer(Num.flwrs ~ Mammal.herb + (1 | Genotype) +
  (1 | Block), data = datExp_Num.flwrs, REML = F)

#Model output
anova(trait.model.6.UnT.Final, type = 3, ddf = "kenward-roger")
rand(trait.model.6.UnT.Final) # Significant Genotype, block and GT:Herb variation


###################################
#### FITNESS ACROSS TREATMENTS ####
###################################

#Note that only final models (post-backward model selection)  are shown

## TOTAL (I.E. ABSOLUTE) FITNESS ##

TF.model.T.Final <- lmer(Total.Seed.mass.T ~ Mammal.herb + HCN +
  Herbivory + Pollination + (1 | Genotype) + (1 | Block) +
  (1 | Herbivory:Genotype) + HCN:Herbivory + HCN:Pollination +
  Herbivory:Pollination + HCN:Herbivory:Pollination, data = datExp,
  REML = F)
summary(TF.model.T.Final)
anova(TF.model.T.Final, type = 3, ddf = "kenward-roger")
rand(TF.model.T.Final)

#Dataset for plotting effects of invertebrate herbivores on absolute fitness
means.TF.Herb <- summarySE(datExp, measurevar = "Total.Seed.mass", groupvars = "Herbivory", na.rm = T)
(3.847784 - 2.604207)/3.847784 # ES

#Dataset for plotting effects of voles on absolute fitness
means.TF.Vole <- summarySE(datExp, measurevar = "Total.Seed.mass", groupvars = "Mammal.herb", na.rm = T)
names(means.TF.Vole)[names(means.TF.Vole) == "Mammal.herb"] <- "Voles"
means.TF.Vole <- within(means.TF.Vole, {
  Voles = ifelse(Voles == 0, "Undamaged", "Damaged")
})
(4.123573 - 1.908990)/4.123573 # ES

#Figure 2C: Absolute fitness by herbivory treatment
means.TF.Herb$Herbivory <- factor(means.TF.Herb$Herbivory, levels=c("Reduced", "Ambient"))
plot.TF.x.Herb <- ggplot(means.TF.Herb,aes(x=Herbivory, y=Total.Seed.mass))+
  geom_errorbar(aes(ymin=Total.Seed.mass-se,ymax=Total.Seed.mass+se),width=0.15,size=0.7)+
  geom_point(size=5.5)+
  xlab("Herbivory")+ylab("Total seed mass (g)") +
  ng1 + theme(axis.text.x=element_text(size=17), axis.text.y=element_text(size=17))

#Figure 2F: Absolute fitness by vole damage.
means.TF.Vole$Voles <- factor(means.TF.Vole$Voles, levels=c("Undamaged", "Damaged"))
plot.TF.x.Vole <- ggplot(means.TF.Vole,aes(x=Voles, y=Total.Seed.mass))+
  geom_errorbar(aes(ymin=Total.Seed.mass-se,ymax=Total.Seed.mass+se),width=0.15,size=0.7)+
  geom_point(size=5.5)+
  xlab("Vole damage")+ylab("Total seed mass (g)") +
  coord_cartesian(ylim = c(1.7, 4.55)) + scale_y_continuous(breaks = seq(from = 2.0, to = 4.5, by = 0.5)) +
  ng1 + theme(axis.text.x=element_text(size=17), axis.text.y=element_text(size=17))

#Save figures 2C and 2F
ggsave("HPS_figures/Figure.2C_TF.x.Herbivory.pdf", plot = plot.TF.x.Herb, width = 5, height = 5, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.2F_TF.x.Voles.pdf", plot = plot.TF.x.Vole, width = 5, height = 5, unit = "in", dpi = 600)

## SEEDS PER INFLORESCENCE ##

PL.model.T.Final <- lmer(Seeds.Inf.T ~ Mammal.herb + HCN + Herbivory +
  Pollination + (1 | Genotype) + (1 | Block) + Mammal.herb:HCN +
  Mammal.herb:Pollination + HCN:Pollination + Mammal.herb:HCN:Pollination,
  data = datExp,
  REML = F)

#Model output
summary(PL.model.T.Final)
anova(PL.model.T.Final, type = 3, ddf = "kenward-roger")
rand(PL.model.T.Final)


###########################################
#### GENOTYPIC SELECTION ANALYSIS: HCN ####
###########################################

## GT SELECTION ANALYSIS LOOKING AT EFFECTS OF TREATMENTS AND HCN ##


#Dataset for genotypic selection analysis with treatments
GTseln.data <- "HPS_data-clean/Santangelo-Thompson-Johnson_Evolution_2018_GTSelnData-all_ExpTreat.txt"
GTSelnData <- read.table(GTseln.data, header = T, fill = T)


#Remove NA's from data
GTSelnData.all <- na.omit(GTSelnData)

#Step 1 -- Global model for genotypic selection analysis (Effects of treatments)
Global.model_GTSeln <- lmer(RF.Seed ~ HCN + Pollination + Herbivory +
  HCN:Pollination + HCN:Herbivory + Pollination:Herbivory + HCN:Herbivory:Pollination +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Biomass.S + Infl.S + Flwrs.S +
  Flwr.date.S:HCN + Flwr.date.S:Pollination + Flwr.date.S:Herbivory +
  Flwr.date.S:HCN:Pollination + Flwr.date.S:HCN:Herbivory +
  Flwr.date.S:Herbivory:Pollination + Flwr.date.S:Herbivory:HCN:Pollination +
  Bnr.wdth.S:HCN + Bnr.wdth.S:Pollination + Bnr.wdth.S:Herbivory +
  Bnr.wdth.S:HCN:Pollination + Bnr.wdth.S:HCN:Herbivory +
  Bnr.wdth.S:Herbivory:Pollination + Bnr.wdth.S:Herbivory:HCN:Pollination +
  Bnr.ht.S:HCN + Bnr.ht.S:Pollination + Bnr.ht.S:Herbivory +
  Bnr.ht.S:HCN:Pollination + Bnr.ht.S:HCN:Herbivory +
  Bnr.ht.S:Herbivory:Pollination + Bnr.ht.S:Herbivory:HCN:Pollination +
  Biomass.S:HCN + Biomass.S:Pollination + Biomass.S:Herbivory +
  Biomass.S:HCN:Pollination + Biomass.S:HCN:Herbivory +
  Biomass.S:Herbivory:Pollination + Biomass.S:Herbivory:HCN:Pollination +
  Infl.S:HCN + Infl.S:Pollination + Infl.S:Herbivory +
  Infl.S:HCN:Pollination + Infl.S:HCN:Herbivory +
  Infl.S:Herbivory:Pollination + Infl.S:Herbivory:HCN:Pollination +
  Flwrs.S:HCN + Flwrs.S:Pollination + Flwrs.S:Herbivory +
  Flwrs.S:HCN:Pollination + Flwrs.S:HCN:Herbivory +
  Flwrs.S:Herbivory:Pollination + Flwrs.S:Herbivory:HCN:Pollination + (1|Genotype),
  data = GTSelnData.all, REML = FALSE)

#Step 2 -- Backward model selection
lmerTest::step(Global.model_GTSeln, ddf = "kenward-roger", type = 3, alpha.random = 0.1, alpha.fixed = 0.05,
              reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, lsmeans.calc = TRUE,
              difflsmeans.calc = TRUE, test.effs = NULL,  keep.effs = "Genotype")

#Step 3 -- Add all main effects back to the model
final.Global.model_GTSeln <- lmer(RF.Seed ~ HCN + Pollination + Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  HCN:Pollination + HCN:Herbivory + HCN:Flwr.date.S + Pollination:Flwr.date.S +
  HCN:Bnr.wdth.S + Herbivory:Bnr.wdth.S + HCN:Bnr.ht.S + Herbivory:Bnr.ht.S +
  HCN:Infl.S + Herbivory:Infl.S + HCN:Pollination:Flwr.date.S +
  HCN:Herbivory:Bnr.wdth.S + HCN:Herbivory:Bnr.ht.S +
  (1 | Genotype), data = GTSelnData.all)

# Selection gradients
summary(final.Global.model_GTSeln)

# P-values
final.Global.model_GTSeln_KR <- anova(final.Global.model_GTSeln, type = 3, ddf = "Kenward-Roger")
final.Global.model_GTSeln_KR

#Write HCN multivariate selection analysis coefficients to dataset and write to disk
final.Global.model_GTSeln_OUT <- broom::tidy(final.Global.model_GTSeln)
final.Global.model_GTSeln_KR_OUT <- broom::tidy(final.Global.model_GTSeln_KR)
write.csv(final.Global.model_GTSeln_OUT, "HPS_tables/Table-1_Gradients_Multivariate-selection_HCN_Treatments.csv")
write.csv(final.Global.model_GTSeln_KR_OUT, "HPS_tables/Table-1_Pvals_Multivariate-selection_HCN_Treatments.csv")

## UNIVARIATE SELECTION GRADIENTS BASED ON SIGNIFICANT INTERACTIONS FROM ABOVE MODEL ##

#Data subsets by HCN
GTSelnData.all.Cyan <- subset(GTSelnData.all, HCN == "Yes")
GTSelnData.all.Acyan <- subset(GTSelnData.all, HCN == "No")

#Reduced model. HCN+ plants
Sel.Cyan <- lmer(RF.Seed ~ Pollination + Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Pollination:Flwr.date.S + Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
  Herbivory:Infl.S + (1 | Genotype),
  data = GTSelnData.all.Cyan)
summary(Sel.Cyan)
anova(Sel.Cyan, type = 3, ddf = "kenward-roger")

#HCN– plants
Sel.Acyan <- lmer(RF.Seed ~ Pollination + Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Pollination:Flwr.date.S + Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
  Herbivory:Infl.S + (1 | Genotype), data = GTSelnData.all.Acyan)
summary(Sel.Acyan)
anova(Sel.Acyan, type = 3, ddf = "kenward-roger")


# Effect size on # inflorescences (Acyan vs. Cyan). 44% weaker among Acyan
(1.39014 - 0.98490)/1.39014

# Number of inflorescences. Fitness residuals. Add to data frame for plotting
GTSelnData.all$InflFitnessResid <- resid(lmer(RF.Seed ~ HCN + Pollination + Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  HCN:Pollination + HCN:Herbivory + HCN:Flwr.date.S + Pollination:Flwr.date.S +
  HCN:Bnr.wdth.S + Herbivory:Bnr.wdth.S + HCN:Bnr.ht.S + Herbivory:Bnr.ht.S +
  HCN:Pollination:Flwr.date.S +
  HCN:Herbivory:Bnr.wdth.S + HCN:Herbivory:Bnr.ht.S +
  (1 | Genotype), data = GTSelnData.all))

#Number of inflorescences residuals. Add to data frame for plotting
GTSelnData.all$InflResid <- resid(lmer(Infl.S ~ HCN + Pollination + Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  HCN:Pollination + HCN:Herbivory + HCN:Flwr.date.S + Pollination:Flwr.date.S +
  HCN:Bnr.wdth.S + Herbivory:Bnr.wdth.S + HCN:Bnr.ht.S + Herbivory:Bnr.ht.S +
  HCN:Pollination:Flwr.date.S +
  HCN:Herbivory:Bnr.wdth.S + HCN:Herbivory:Bnr.ht.S +
  (1 | Genotype), data = GTSelnData.all))

#Figure 3D. Plot number of inflorescences by HCN. 3 points not shown
GTSelnData.all[order(GTSelnData.all$InflFitnessResid,decreasing=T)[1:5],]
PlotInfl.x.HCN <- ggplot(GTSelnData.all, aes(x = InflResid, y = InflFitnessResid, group = HCN)) +
  labs(x = "Number of inflorescences (residual)", y = "Relative fitness (residual)") +
  geom_point(aes(shape = HCN, fill = HCN)) +
  scale_shape_manual(labels = c("HCN-","HCN+"), values=c(24, 21)) +
  scale_fill_manual(labels = c("HCN-","HCN+"), values=c("white", "black")) +
  geom_smooth(method = "lm", se = F, colour = "black", size = 1.05, fullrange = T, aes(linetype = HCN))  +
  scale_linetype_manual(labels = c("HCN-","HCN+"), values=c(2, 1)) +
  coord_cartesian(ylim = c(-2.3, 2.0)) + scale_y_continuous(breaks = seq(from = -2, to = 2.0, by = 0.5)) +
  ng1 + theme(legend.title = element_blank())

#Save figures 3D
ggsave("HPS_figures/Figure.3D_Sel.Infl.x.HCN.pdf", plot = PlotInfl.x.HCN, width = 5, height = 5, unit = "in", dpi = 600)

#Data subset by herbivory treatment
GTSelnData.all.Am <- subset(GTSelnData.all, Herbivory == "Ambient")
GTSelnData.all.Red <- subset(GTSelnData.all, Herbivory == "Reduced")

#Ambient herbivory
Sel.Am <- lmer(RF.Seed ~ HCN + Pollination +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  HCN:Pollination + HCN:Flwr.date.S + Pollination:Flwr.date.S +
  HCN:Bnr.wdth.S + HCN:Bnr.ht.S + HCN:Infl.S + HCN:Pollination:Flwr.date.S +
  (1 | Genotype), data = GTSelnData.all.Am)
summary(Sel.Am)
anova(Sel.Am, type = 3, ddf = "kenward-roger")

#Reduced herbivory
Sel.Red <- lmer(RF.Seed ~ HCN + Pollination +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  HCN:Pollination + HCN:Flwr.date.S + Pollination:Flwr.date.S +
  HCN:Bnr.wdth.S + HCN:Bnr.ht.S + HCN:Infl.S + HCN:Pollination:Flwr.date.S +
  (1 | Genotype), data = GTSelnData.all.Red)
summary(Sel.Red)
anova(Sel.Red, type = 3, ddf = "kenward-roger")

# Effect size on # of inflorescences (Ambient vs. Reduced). 28% weaker for Ambient
(1.343408 - 0.914871)/1.343408

#Figure 3A. Plot number of inflorescences by Herbivory. 2 points not shown
PlotInfl.x.Herb <- ggplot(GTSelnData.all, aes(x = InflResid, y = InflFitnessResid, group = Herbivory)) +
  labs(x = "Number of inflorescences (residual)", y = "Relative fitness (residual)") +
  geom_point(aes(shape = Herbivory, fill = Herbivory)) +
  scale_shape_manual(labels = c("Ambient herbivory","Reduced herbivory"), values=c(22, 23)) +
  scale_fill_manual(labels = c("Ambient herbivory","Reduced herbivory"), values=c("white", "black")) +
  geom_smooth(method = "lm", se = F, colour = "black", size = 1.05, fullrange = T, aes(linetype = Herbivory))  +
  scale_linetype_manual(labels = c("Ambient herbivory","Reduced herbivory"), values=c(2, 1)) +
  coord_cartesian(ylim = c(-2.3, 2.6)) + scale_y_continuous(breaks = seq(from = -2, to = 2.5, by = 0.5)) +
  ng1 + theme(legend.title = element_blank())

ggsave("HPS_figures/Figure.3A_Sel.Infl.x.Herb.pdf", plot = PlotInfl.x.Herb, width = 5, height = 5, unit = "in", dpi = 600)

## 3-WAY INTERACTION BETWEEN HCN, HERBIVORY TREATMEN AND BANNER LENGTH/WIDTH

#Data subsets by HCN and herbivory treatment
GTSelnData.all.Cyan.Am <- subset(GTSelnData.all, HCN == "Yes" & Herbivory =="Ambient")
GTSelnData.all.Acyan.Am <- subset(GTSelnData.all, HCN == "No" & Herbivory =="Ambient")
GTSelnData.all.Cyan.Red <- subset(GTSelnData.all, HCN == "Yes" & Herbivory =="Reduced")
GTSelnData.all.Acyan.Red <- subset(GTSelnData.all, HCN == "No" & Herbivory =="Reduced")

#HCN+, Ambient herbivory
top.model.Cyan.Am <- lmer(RF.Seed ~ Pollination +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Pollination:Flwr.date.S +
  (1 | Genotype), data = GTSelnData.all.Cyan.Am)
summary(top.model.Cyan.Am)
anova(top.model.Cyan.Am, type = 3, ddf = "kenward-roger")

#HCN–, Ambient herbivory
top.model.Acyan.Am <- lmer(RF.Seed ~ Pollination +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Pollination:Flwr.date.S +
  (1 | Genotype), data = GTSelnData.all.Acyan.Am)
summary(top.model.Acyan.Am)
anova(top.model.Acyan.Am, type = 3, ddf = "kenward-roger")

#HCN+, Reduced herbivory
top.model.Cyan.Red <- lmer(RF.Seed ~ Pollination +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Pollination:Flwr.date.S +
  (1 | Genotype), data = GTSelnData.all.Cyan.Red)
summary(top.model.Cyan.Red)
anova(top.model.Cyan.Red, type = 3, ddf = "kenward-roger")

#HCN–, Reduced herbivory
top.model.Acyan.Red <- lmer(RF.Seed ~ Pollination +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Pollination:Flwr.date.S +
  (1 | Genotype), data = GTSelnData.all.Acyan.Red)
summary(top.model.Acyan.Red)
anova(top.model.Acyan.Red, type = 3, ddf = "kenward-roger")

#Dataset for Figure 4A in paper.
Coefs <- rbind(coef(summary(top.model.Cyan.Am))[c("Bnr.wdth.S", "Bnr.ht.S"),1:2],
              coef(summary(top.model.Acyan.Am))[c("Bnr.wdth.S", "Bnr.ht.S"),1:2],
              coef(summary(top.model.Cyan.Red))[c("Bnr.wdth.S", "Bnr.ht.S"),1:2],
              coef(summary(top.model.Acyan.Red))[c("Bnr.wdth.S", "Bnr.ht.S"),1:2])
Herbivory <-  c("Ambient", "Ambient", "Ambient", "Ambient", "Reduced", "Reduced", "Reduced", "Reduced")
HCN <- c("Yes", "Yes", "No", "No", "Yes", "Yes", "No", "No")
Bnr.Size.Gradients <- as.data.frame(cbind(Coefs, Herbivory, HCN))
names(Bnr.Size.Gradients)[names(Bnr.Size.Gradients) == "Std. Error"] <- "se"
Bnr.Size.Gradients$Estimate <- as.numeric(as.character(Bnr.Size.Gradients$Estimate))
Bnr.Size.Gradients$se <- as.numeric(as.character(Bnr.Size.Gradients$se))
Bnr.Size.Gradients$Trait <- rownames(Bnr.Size.Gradients)
Bnr.Size.Gradients$ci <- 1.96*Bnr.Size.Gradients$se

#Figure 4A in paper. Selection on banner width by herbivory and HCN
Bnr.wdth.Gradients <- subset(Bnr.Size.Gradients, Trait == "Bnr.wdth.S")
plotBnr.wdth.Gradients <- ggplot(Bnr.wdth.Gradients, aes(x = HCN, y = Estimate, shape = Herbivory, fill = Herbivory)) +
  labs(x = "HCN", y = "Selection gradient") +
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=Estimate-ci,ymax=Estimate+ci),width=0.15,size=0.7,position = position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5), size = 4.5, color = "black") +
  scale_shape_manual(labels = c("Ambient", "Reduced"),values=c(22, 23)) +
  scale_fill_manual(labels = c("Ambient", "Reduced"),values=c("white", "black")) +
  coord_cartesian(ylim = c(-0.45, 0.40)) + scale_y_continuous(breaks = seq(from = -0.40, to = 0.40, by = 0.1)) +
  ng1 + theme(axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())

#Figure 4B in paper. Selection on banner length by herbivory and HCN
Bnr.lgth.Gradients <- subset(Bnr.Size.Gradients, Trait == "Bnr.ht.S")
scaleFUN <- function(x) sprintf("%.1f", x)
plotBnr.lgth.Gradients <- ggplot(Bnr.lgth.Gradients, aes(x = HCN, y = Estimate, shape = Herbivory, fill = Herbivory)) +
  labs(x = "HCN", y = "Selection gradient") +
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=Estimate-ci,ymax=Estimate+ci),width=0.15,size=0.7,position = position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5), size = 4.5, color = "black") +
  scale_shape_manual(labels = c("Ambient", "Reduced"),values=c(22, 23)) +
  scale_fill_manual(labels = c("Ambient", "Reduced"),values=c("white", "black")) +
  coord_cartesian(ylim = c(-0.6, 0.40)) + scale_y_continuous(breaks = seq(from = -0.60, to = 0.40, by = 0.1), labels = scaleFUN) +
  ng1 + theme(axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())

#Save figures 4A and 4B
ggsave("HPS_figures/Figure.4A_Sel.BW.x.HCN.x.Herb.pdf", plot = plotBnr.wdth.Gradients, width = 5, height = 5, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.4B_Sel.BL.x.HCN.x.Herb.pdf", plot = plotBnr.lgth.Gradients, width = 5, height = 5, unit = "in", dpi = 600)

## 3-WAY INTERACTION BETWEEN HCN, POLLINATION TREATMENT AND DATE TO FIRST FLOWER

#Data subsets by HCN and pollination treatment
GTSelnData.all.Cyan.Open <- subset(GTSelnData.all, HCN == "Yes" & Pollination =="Open")
GTSelnData.all.Acyan.Open <- subset(GTSelnData.all, HCN == "No" & Pollination =="Open")
GTSelnData.all.Cyan.Supp <- subset(GTSelnData.all, HCN == "Yes" & Pollination =="Supp")
GTSelnData.all.Acyan.Supp <- subset(GTSelnData.all, HCN == "No" & Pollination =="Supp")

#HCN+, Open pollination
top.model.Cyan.Open <- lmer(RF.Seed ~ Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
  Herbivory:Infl.S + (1 | Genotype),
  data = GTSelnData.all.Cyan.Open)
summary(top.model.Cyan.Open)
anova(top.model.Cyan.Open, type = 3, ddf = "kenward-roger")

#HCN–, Open pollination
top.model.Acyan.Open <- lmer(RF.Seed ~ Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
  Herbivory:Infl.S + (1 | Genotype),
  data = GTSelnData.all.Acyan.Open)
summary(top.model.Acyan.Open)
anova(top.model.Acyan.Open, type = 3, ddf = "kenward-roger")


#HCN+, Supplemental pollination
top.model.Cyan.Supp <- lmer(RF.Seed ~ Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
  Herbivory:Infl.S + (1 | Genotype),
  data = GTSelnData.all.Cyan.Supp)
summary(top.model.Cyan.Supp)
anova(top.model.Cyan.Supp, type = 3, ddf = "kenward-roger")

#HCN–, Supplemental pollination
top.model.Acyan.Supp <- lmer(RF.Seed ~ Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
  Herbivory:Infl.S + (1 | Genotype),
  data = GTSelnData.all.Acyan.Supp)
summary(top.model.Acyan.Supp)
anova(top.model.Acyan.Supp, type = 3, ddf = "kenward-roger")

#Dataset for Figure S10 in paper.
Coefs <- rbind(coef(summary(top.model.Cyan.Open))["Flwr.date.S",1:2],
              coef(summary(top.model.Acyan.Open))["Flwr.date.S",1:2],
              coef(summary(top.model.Cyan.Supp))["Flwr.date.S",1:2],
              coef(summary(top.model.Acyan.Supp))["Flwr.date.S",1:2])
Pollination <-  c("Open", "Open", "Supplemental", "Supplemental")
HCN <- c("Yes", "No", "Yes", "No")
Flwr.date.Gradients <- as.data.frame(cbind(Coefs, Pollination, HCN))
names(Flwr.date.Gradients)[names(Flwr.date.Gradients) == "Std. Error"] <- "se"
Flwr.date.Gradients$Estimate <- as.numeric(as.character(Flwr.date.Gradients$Estimate))
Flwr.date.Gradients$se <- as.numeric(as.character(Flwr.date.Gradients$se))
Flwr.date.Gradients$Trait <- "Date to first flower"
Flwr.date.Gradients$ci <- 1.96*Flwr.date.Gradients$se

#Figure S10 in paper. Selection on banner width by herbivory and HCN
Flwr.date.Gradients$Pollination <- factor(Flwr.date.Gradients$Pollination, levels=c("Open", "Supplemental"))
Flwr.date.Gradients$HCN <- factor(Flwr.date.Gradients$HCN, levels=c("No", "Yes"))
plotFlwr.date.Gradients <- ggplot(Flwr.date.Gradients, aes(x = HCN, y = Estimate, group = Pollination, shape = Pollination, fill = Pollination)) +
  labs(x = "HCN", y = "Selection gradient") +
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=Estimate-ci,ymax=Estimate+ci),width=0.15,size=0.7,position = position_dodge(width = 0.5)) +
  geom_point(position=position_dodge(width = 0.5), size = 4.5, color = "black") +
  scale_shape_manual(labels = c("Open", "Supplemental"),values=c(22, 23)) +
  scale_fill_manual(labels = c("Open", "Supplemental"),values=c("white", "black")) +
  # coord_cartesian(ylim = c(-0.45, 0.40)) + scale_y_continuous(breaks = seq(from = -0.40, to = 0.40, by = 0.1)) +
  ng1 + theme(axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())

ggsave("HPS_figures/Figure.S10_Sel.FF.x.HCN.x.Poll.pdf", plot = plotFlwr.date.Gradients, width = 5, height = 5, unit = "in", dpi = 600)


#############################################
#### GENOTYPIC SELECTION ANALYSIS: VOLES ####
#############################################

#Dataset for genotypic selection analysis of vole damage (name = GTSelnData-all_Voles.txt)
GTseln.Voles <- "HPS_data-clean/Santangelo-Thompson-Johnson_Evolution_2018_GTSelnData-all_Voles.txt"
GTSelnData.all.Voles <- read.table(GTseln.Voles, header = T, fill = T)

GTSelnData.all.Voles <- na.omit(GTSelnData.all.Voles)

#Global model for genotypic selection analysis (Effects of voles)
Global.model_GTSeln.Voles <- lmer(RF.Seed ~ Mammal.herb  + HCN + Mammal.herb:HCN +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  Biomass.S + Infl.S + Flwrs.S +
  HCN:Flwr.date.S + HCN:Bnr.ht.S + HCN:Bnr.wdth.S +
  HCN:Biomass.S + HCN:Infl.S + HCN:Flwrs.S +
  Mammal.herb:Biomass.S + Mammal.herb:Bnr.ht.S +
  Mammal.herb:Bnr.ht.S + Mammal.herb:Bnr.wdth.S +
  Mammal.herb:Infl.S + Mammal.herb:Flwrs.S +
  HCN:Mammal.herb:Biomass.S + HCN:Mammal.herb:Bnr.ht.S +
  HCN:Mammal.herb:Bnr.ht.S + HCN:Mammal.herb:Bnr.wdth.S +
  HCN:Mammal.herb:Infl.S + HCN:Mammal.herb:Flwrs.S +
    (1|Genotype),
  data = GTSelnData.all.Voles, REML = FALSE)

#Backward model selection of above model
lmerTest::step(Global.model_GTSeln.Voles, ddf = "kenward-roger", type = 3, alpha.random = 0.1, alpha.fixed = 0.05,
              reduce.fixed = TRUE, reduce.random = FALSE, fixed.calc = TRUE, lsmeans.calc = TRUE,
              difflsmeans.calc = TRUE, test.effs = NULL,  keep.effs = "Genotype")

#Optimized model from backward model selection
final.model_GTSeln.Voles <- lmer(RF.Seed ~ Mammal.herb + HCN + Bnr.wdth.S +
   Infl.S + Flwrs.S + Mammal.herb:HCN + HCN:Bnr.wdth.S +
   HCN:Infl.S + Mammal.herb:Bnr.wdth.S + Mammal.herb:Infl.S +
   Mammal.herb:Flwrs.S + Mammal.herb:HCN:Bnr.wdth.S +
   (1 | Genotype), data = GTSelnData.all.Voles)

# Selection gradients
summary(final.model_GTSeln.Voles)

# P-values
final.model_GTSeln.Voles_KR <- anova(final.model_GTSeln.Voles, type = 3, ddf = "Kenward-Roger")
final.model_GTSeln.Voles_KR

#Write HCN multivariate selection analysis coefficients to dataset
final.model_GTSeln.Voles_OUT <- broom::tidy(final.model_GTSeln.Voles)
final.model_GTSeln.Voles_KR_OUT <- broom::tidy(final.model_GTSeln.Voles_KR)
write.csv(final.model_GTSeln.Voles_OUT, "HPS_tables/Table-S5_Gradients_Multivariate-selection_Voles_Treatments.csv")
write.csv(final.model_GTSeln.Voles_KR_OUT, "HPS_tables/Table-S5_Pvals_Multivariate-selection_Voles_Treatments.csv")

## UNIVARIATE SELECTION GRADIENTS FROM SIGNIFICANT INTERACTION IN GT SELECTION ANALYSIS OF VOLES ##

#Data subsets by vole damage
GTSelnData.all.Voles.Dmg <- subset(GTSelnData.all.Voles, Mammal.herb == "1")
GTSelnData.all.Voles.Undmg <- subset(GTSelnData.all.Voles, Mammal.herb == "0")

#Damaged by voles. Genotype not included as random effect since only one obs. per genotype
top.model.Voles.Dmg <- lm(RF.Seed ~ HCN + Bnr.wdth.S +
  Infl.S + Flwrs.S + HCN:Bnr.wdth.S +
  HCN:Infl.S, data = GTSelnData.all.Voles.Dmg)
summary(top.model.Voles.Dmg)
car::Anova(top.model.Voles.Dmg, type = 3, ddf = "kenward-roger")

#Undamaged by voles
top.model.Voles.Undmg <- lm(RF.Seed ~ HCN + Bnr.wdth.S +
  Infl.S + Flwrs.S + HCN:Bnr.wdth.S +
  HCN:Infl.S, data = GTSelnData.all.Voles.Undmg)
summary(top.model.Voles.Undmg)
car::Anova(top.model.Voles.Undmg, type = 3, ddf = "kenward-roger")

# Effect size on # of Inflorescences (Damaged vs. Undamaged). 94% weaker.
(1.21121 - 0.67866)/1.21121
# Effect size on # of flowers (Damaged vs. Undamaged). 264% weaker.
(0.23387 - 0.06416)/0.23387

#Data subsets by vole damage and HCN
GTSelnData.all.Voles.Cyan.Dmg <- subset(GTSelnData.all.Voles, HCN == "Yes" & Mammal.herb == "1")
GTSelnData.all.Voles.Acyan.Dmg <- subset(GTSelnData.all.Voles, HCN == "No" & Mammal.herb == "1")
GTSelnData.all.Voles.Cyan.Undmg <- subset(GTSelnData.all.Voles, HCN == "Yes" & Mammal.herb == "0")
GTSelnData.all.Voles.Acyan.Undmg <- subset(GTSelnData.all.Voles, HCN == "No" & Mammal.herb == "0")

# HCN+, Damaged
top.model.Cyan.Dmg <- lm(RF.Seed ~ Bnr.wdth.S +
  Infl.S + Flwrs.S,
  data = GTSelnData.all.Voles)
summary(top.model.Cyan.Dmg)

# HCN-, Damaged
top.model.Acyan.Dmg <- lm(RF.Seed ~ Bnr.wdth.S +
  Infl.S + Flwrs.S,
  data = GTSelnData.all.Voles.Acyan.Dmg)
summary(top.model.Acyan.Dmg)

# HCN+, Undamaged
top.model.Cyan.Dmg <- lm(RF.Seed ~ Bnr.wdth.S +
  Infl.S + Flwrs.S,
  data = GTSelnData.all.Voles.Cyan.Undmg)
summary(top.model.Cyan.Dmg)

# HCN-, Undamaged
top.model.Acyan.Undmg <- lm(RF.Seed ~ Bnr.wdth.S +
  Infl.S + Flwrs.S,
  data = GTSelnData.all.Voles.Acyan.Undmg)
summary(top.model.Acyan.Undmg)

#Number of Inflorescences. Fitness residuals. Add to dataframe
GTSelnData.all.Voles$InflFitnessResid <- resid(lmer(RF.Seed ~ Mammal.herb + HCN + Bnr.wdth.S +
  Flwrs.S + Mammal.herb:HCN + HCN:Bnr.wdth.S +
  Mammal.herb:Bnr.wdth.S +
  Mammal.herb:Flwrs.S + Mammal.herb:HCN:Bnr.wdth.S +
  (1 | Genotype), data = GTSelnData.all.Voles))

#Number of inflorescences residuals. Add to dataframe
GTSelnData.all.Voles$InflResid <- resid(lmer(Infl.S ~ Mammal.herb + HCN + Bnr.wdth.S +
  Flwrs.S + Mammal.herb:HCN + HCN:Bnr.wdth.S +
  Mammal.herb:Bnr.wdth.S +
  Mammal.herb:Flwrs.S + Mammal.herb:HCN:Bnr.wdth.S +
  (1 | Genotype), data = GTSelnData.all.Voles))

#Figure 3B. Plot number of inflorescences by vole damage.
GTSelnData.all.Voles[order(GTSelnData.all.Voles$InflFitnessResid,decreasing=T)[1:5],]
PlotInfl.x.Voles <- ggplot(GTSelnData.all.Voles, aes(x = InflResid, y = InflFitnessResid, group = factor(Mammal.herb))) +
  labs(x = "Number of inflorescences (residual)", y = "Relative fitness (residual)") +
  geom_point(aes(shape = factor(Mammal.herb), fill = factor(Mammal.herb))) +
  scale_shape_manual(labels = c("No vole damage","Vole damage"), values=c(21, 25)) +
  scale_fill_manual(labels = c("No vole damage","Vole damage"), values=c("black", "white")) +
  geom_smooth(method = "lm", se = F, colour = "black", size = 1.05, fullrange = T, aes(linetype = factor(Mammal.herb)))  +
  scale_linetype_manual(labels = c("No vole damage","Vole damage"), values=c(1, 2)) +
  coord_cartesian(ylim = c(-1.5, 2.0)) + scale_y_continuous(breaks = seq(from = -2, to = 3.0, by = 0.5)) +
  ng1 + theme(legend.title = element_blank())

ggsave("HPS_figures/Figure.3B_Sel.Infl.x.Voles.pdf", plot = PlotInfl.x.Voles, width = 5, height = 5, unit = "in", dpi = 600)

#Number of flowers Fitness residuals. Add to dataframe
GTSelnData.all.Voles$FlwrsFitnessResid <- resid(lmer(RF.Seed ~ Mammal.herb + HCN + Bnr.wdth.S +
  Infl.S + Mammal.herb:HCN + HCN:Bnr.wdth.S +
  HCN:Infl.S + Mammal.herb:Bnr.wdth.S + Mammal.herb:Infl.S +
  Mammal.herb:HCN:Bnr.wdth.S +
  (1 | Genotype), data = GTSelnData.all.Voles))

#Number of flowers residuals. Add to dataframe
GTSelnData.all.Voles$FlwrsResid <- resid(lmer(Flwrs.S ~ Mammal.herb + HCN + Bnr.wdth.S +
  Infl.S + Mammal.herb:HCN + HCN:Bnr.wdth.S +
  HCN:Infl.S + Mammal.herb:Bnr.wdth.S + Mammal.herb:Infl.S +
  Mammal.herb:HCN:Bnr.wdth.S +
  (1 | Genotype), data = GTSelnData.all.Voles))

#Figure 3C. Plot number of flowers per inflorescence by vole damage. 3 point not shown
GTSelnData.all.Voles[order(GTSelnData.all.Voles$FlwrsFitnessResid,decreasing=T)[1:5],]
PlotFlwrs.x.Voles <- ggplot(GTSelnData.all.Voles, aes(x = FlwrsResid, y = FlwrsFitnessResid, group = factor(Mammal.herb))) +
  labs(x = "Number of flowers (residual)", y = "Relative fitness (residual)") +
  geom_point(aes(shape = factor(Mammal.herb), fill = factor(Mammal.herb))) +
  scale_shape_manual(labels = c("No vole damage","Vole damage"), values=c(21, 25)) +
  scale_fill_manual(labels = c("No vole damage","Vole damage"), values=c("black", "white")) +
  geom_smooth(method = "lm", se = F, colour = "black", size = 1.05, fullrange = T, aes(linetype = factor(Mammal.herb)))  +
  scale_linetype_manual(labels = c("No vole damage","Vole damage"), values=c(1, 2)) +
  coord_cartesian(ylim = c(-0.55, 0.5)) + scale_y_continuous(breaks = seq(from = -0.55, to = 0.5, by = 0.2)) +
  ng1 + theme(legend.title = element_blank())

ggsave("HPS_figures/Figure.3C_Sel.Flwrs.x.Voles.pdf", plot = PlotFlwrs.x.Voles, width = 5, height = 5, unit = "in", dpi = 600)



###########################################
#### GENOTYPIC SELECTION ANALYSIS: Ac ####
###########################################

# Here we see if altered patterns of selection due to the expression of defence
# is due to allocation costs or selection influenced by investment in either
# component gene.

## GT SELECTION ANALYSIS LOOKING AT EFFECTS OF TREATMENTS ##

#Dataset for genotypic selection analysis with treatments
GTSelnData.all.Acy <- subset(GTSelnData.all, HCN == "No")

#Step 1 -- Reduced from multivariate selection analysis of HCN
final.Global.model_GTSeln_Ac <- lmer(RF.Seed ~ Glycosides.Ac + Pollination + Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Glycosides.Ac:Pollination + Glycosides.Ac:Herbivory + Glycosides.Ac:Flwr.date.S + Pollination:Flwr.date.S +
  Glycosides.Ac:Bnr.wdth.S + Herbivory:Bnr.wdth.S + Glycosides.Ac:Bnr.ht.S + Herbivory:Bnr.ht.S +
  Glycosides.Ac:Infl.S + Herbivory:Infl.S + Glycosides.Ac:Pollination:Flwr.date.S +
  Glycosides.Ac:Herbivory:Bnr.wdth.S + Glycosides.Ac:Herbivory:Bnr.ht.S +
  (1 | Genotype), data = GTSelnData.all.Acy)

# Selection gradients
summary(final.Global.model_GTSeln_Ac)

# P-values
final.Global.model_GTSeln_Ac_KR <- anova(final.Global.model_GTSeln_Ac, type = 3, ddf = "Kenward-Roger")
final.Global.model_GTSeln_Ac_KR

#Write HCN multivariate selection analysis coefficients to dataset
final.Global.model_GTSeln_Ac_OUT <- broom::tidy(final.Global.model_GTSeln_Ac)
final.Global.model_GTSeln_Ac_KR_OUT <- broom::tidy(final.Global.model_GTSeln_Ac_KR)
write.csv(final.Global.model_GTSeln_Ac_OUT, "HPS_tables/Table-S6_Gradients_Multivariate-selection_Ac_Treatments.csv")
write.csv(final.Global.model_GTSeln_Ac_KR_OUT, "HPS_tables/Table-S6_Pvals_Multivariate-selection_Ac_Treatments.csv")

## 3-WAY INTERACTION BETWEEN CYP79D15, POLLINATION AND DATE TO FIRST FLOWER

#Data subsets by HCN and pollination treatment
GTSelnData.all.CYP.Open <- subset(GTSelnData.all.Acy, Glycosides.Ac == "Yes" & Pollination =="Open")
GTSelnData.all.ACYP.Open <- subset(GTSelnData.all.Acy, Glycosides.Ac == "No" & Pollination =="Open")
GTSelnData.all.CYP.Supp <- subset(GTSelnData.all.Acy, Glycosides.Ac == "Yes" & Pollination =="Supp")
GTSelnData.all.ACYP.Supp <- subset(GTSelnData.all.Acy, Glycosides.Ac == "No" & Pollination =="Supp")

#CYP+, Open pollination
top.model.CYP.Open <- lmer(RF.Seed ~ Herbivory +
                             Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
                             Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
                             Herbivory:Infl.S + (1 | Genotype),
                            data = GTSelnData.all.CYP.Open)
summary(top.model.CYP.Open)
anova(top.model.CYP.Open, type = 3, ddf = "kenward-roger")

#CYP–, Open pollination
top.model.ACYP.Open <- lmer(RF.Seed ~ Herbivory +
                               Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
                               Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
                               Herbivory:Infl.S + (1 | Genotype),
                             data = GTSelnData.all.ACYP.Open)
summary(top.model.ACYP.Open)
anova(top.model.ACYP.Open, type = 3, ddf = "kenward-roger")


#CYP+, Supplemental pollination
top.model.CYP.Supp <- lmer(RF.Seed ~ Herbivory +
                              Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
                              Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
                              Herbivory:Infl.S + (1 | Genotype),
                            data = GTSelnData.all.CYP.Supp)
summary(top.model.CYP.Supp)
anova(top.model.CYP.Supp, type = 3, ddf = "kenward-roger")

#CYP–, Supplemental pollination
top.model.ACYP.Supp <- lmer(RF.Seed ~ Herbivory +
                               Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
                               Herbivory:Bnr.wdth.S + Herbivory:Bnr.ht.S +
                               Herbivory:Infl.S + (1 | Genotype),
                             data = GTSelnData.all.ACYP.Supp)
summary(top.model.ACYP.Supp)
anova(top.model.ACYP.Supp, type = 3, ddf = "kenward-roger")

###########################################
#### GENOTYPIC SELECTION ANALYSIS: Li ####
###########################################

## GT SELECTION ANALYSIS LOOKING AT EFFECTS OF TREATMENTS ##

#Step 1 -- Reduced from multivariate selection analysis of HCN
final.Global.model_GTSeln_Li <- lmer(RF.Seed ~ Linamarase.Li + Pollination + Herbivory +
  Flwr.date.S + Bnr.ht.S + Bnr.wdth.S + Infl.S +
  Linamarase.Li:Pollination + Linamarase.Li:Herbivory + Linamarase.Li:Flwr.date.S + Pollination:Flwr.date.S +
  Linamarase.Li:Bnr.wdth.S + Herbivory:Bnr.wdth.S + Linamarase.Li:Bnr.ht.S + Herbivory:Bnr.ht.S +
  Linamarase.Li:Infl.S + Herbivory:Infl.S + Linamarase.Li:Pollination:Flwr.date.S +
  Linamarase.Li:Herbivory:Bnr.wdth.S + Linamarase.Li:Herbivory:Bnr.ht.S +
  (1 | Genotype), data = GTSelnData.all.Acy)

# Selection gradients
summary(final.Global.model_GTSeln_Li)

# P-values
final.Global.model_GTSeln_Li_KR <- anova(final.Global.model_GTSeln_Li, type = 3, ddf = "Kenward-Roger")
final.Global.model_GTSeln_Li_KR

#Write HCN multivariate selection analysis coefficients to dataset
final.Global.model_GTSeln_Li_OUT <- broom::tidy(final.Global.model_GTSeln_Li)
final.Global.model_GTSeln_Li_KR_OUT <- broom::tidy(final.Global.model_GTSeln_Li_KR)
write.csv(final.Global.model_GTSeln_Li_OUT, "HPS_tables/Table-S7_Gradients_Multivariate-selection_Li_Treatments.csv")
write.csv(final.Global.model_GTSeln_Li_KR_OUT, "HPS_tables/Table-S7_Pvals_Multivariate-selection_Li_Treatments.csv")

##########################################################
#### FIGURE SUMMARIZING SELECTION BY DIFFERENT AGENTS ####
##########################################################

## AMBIENT HERBIVORY ##
GTSelnData.all.HerbAmb <- subset(GTSelnData.all, Herbivory == "Ambient")
Global.model_GTSeln.HerbAmb <- lmer(RF.Seed ~ Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  Biomass.S + Infl.S + Flwrs.S +
  (1|Genotype), data = GTSelnData.all.HerbAmb)
Coeffs.HerbAmb <- as.data.frame(coef(summary(Global.model_GTSeln.HerbAmb))[ , "Estimate"])
names(Coeffs.HerbAmb)[names(Coeffs.HerbAmb) == "coef(summary(Global.model_GTSeln.HerbAmb))[, \"Estimate\"]"] <- "Gradient.HerbAmb"
Coeffs.HerbAmb$Trait <- rownames(Coeffs.HerbAmb)

## REDUCED HERBIVORY ##
GTSelnData.all.HerbRed <- subset(GTSelnData.all, Herbivory == "Reduced")
Global.model_GTSeln.HerbRed <- lmer(RF.Seed ~ Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  Biomass.S + Infl.S + Flwrs.S +
  (1|Genotype), data = GTSelnData.all.HerbRed)
Coeffs.HerbRed <- as.data.frame(coef(summary(Global.model_GTSeln.HerbRed))[ , "Estimate"])
names(Coeffs.HerbRed)[names(Coeffs.HerbRed) == "coef(summary(Global.model_GTSeln.HerbRed))[, \"Estimate\"]"] <- "Gradient.HerbRed"
Coeffs.HerbRed$Trait <- rownames(Coeffs.HerbRed)

## OPEN POLLINATION ##
GTSelnData.all.PollOp <- subset(GTSelnData.all, Pollination == "Open")
Global.model_GTSeln.PollOp <- lmer(RF.Seed ~ Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  Biomass.S + Infl.S + Flwrs.S +
  (1|Genotype), data = GTSelnData.all.PollOp)
Coeffs.PollOp <- as.data.frame(coef(summary(Global.model_GTSeln.PollOp))[ , "Estimate"])
names(Coeffs.PollOp)[names(Coeffs.PollOp) == "coef(summary(Global.model_GTSeln.PollOp))[, \"Estimate\"]"] <- "Gradient.PollOp"
Coeffs.PollOp$Trait <- rownames(Coeffs.PollOp)

## SUPPLEMENTAL POLLINATION ##
GTSelnData.all.PollSupp <- subset(GTSelnData.all, Pollination == "Supp")
Global.model_GTSeln.PollSupp <- lmer(RF.Seed ~ Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  Biomass.S + Infl.S + Flwrs.S +
  (1|Genotype), data = GTSelnData.all.PollSupp)
Coeffs.PollSupp <- as.data.frame(coef(summary(Global.model_GTSeln.PollSupp))[ , "Estimate"])
names(Coeffs.PollSupp)[names(Coeffs.PollSupp) == "coef(summary(Global.model_GTSeln.PollSupp))[, \"Estimate\"]"] <- "Gradient.PollSupp"
Coeffs.PollSupp$Trait <- rownames(Coeffs.PollSupp)

## HCN+ ##
GTSelnData.all.Cyan <- subset(GTSelnData.all, HCN == "Yes")
Global.model_GTSeln.Cyan <- lmer(RF.Seed ~ Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  Biomass.S + Infl.S + Flwrs.S +
  (1|Genotype), data = GTSelnData.all.Cyan)
Coeffs.Cyan <- as.data.frame(coef(summary(Global.model_GTSeln.Cyan))[ , "Estimate"])
names(Coeffs.Cyan)[names(Coeffs.Cyan) == "coef(summary(Global.model_GTSeln.Cyan))[, \"Estimate\"]"] <- "Gradient.Cyan"
Coeffs.Cyan$Trait <- rownames(Coeffs.Cyan)

## HCN- ##
GTSelnData.all.Acyan <- subset(GTSelnData.all, HCN == "No")
Global.model_GTSeln.Acyan <- lmer(RF.Seed ~ Flwr.date.S + Bnr.ht.S + Bnr.wdth.S +
  Biomass.S + Infl.S + Flwrs.S +
  (1|Genotype), data = GTSelnData.all.Acyan)
Coeffs.Acyan <- as.data.frame(coef(summary(Global.model_GTSeln.Acyan))[ , "Estimate"])
names(Coeffs.Acyan)[names(Coeffs.Acyan) == "coef(summary(Global.model_GTSeln.Acyan))[, \"Estimate\"]"] <- "Gradient.Acyan"
Coeffs.Acyan$Trait <- rownames(Coeffs.Acyan)

Sel.Grad <- Reduce(function(...) merge(..., by = "Trait"),
                   list(Coeffs.HerbAmb, Coeffs.HerbRed, Coeffs.PollOp, Coeffs.PollSupp, Coeffs.Cyan, Coeffs.Acyan))
Sel.Grad <- Sel.Grad[-which(Sel.Grad$Trait == "(Intercept)"),]

Sel.Grad$HerbMedSel <- Sel.Grad$Gradient.HerbAmb - Sel.Grad$Gradient.HerbRed
Sel.Grad$PollMedSel <- Sel.Grad$Gradient.PollOp - Sel.Grad$Gradient.PollSupp
Sel.Grad$DefMedSel <- Sel.Grad$Gradient.Cyan - Sel.Grad$Gradient.Acyan


AgentMedSel <- dplyr::select(Sel.Grad, Trait, HerbMedSel, PollMedSel, DefMedSel)
AgentMedSel <- data.table::melt(AgentMedSel,id.vars = c("Trait"),
                                variable.name = "Agent",
                                value.name = "Gradient")
AgentMedSel <- within(AgentMedSel, {
  Agent = ifelse(Agent == "HerbMedSel","Herbivore",
                 ifelse(Agent == "PollMedSel", "Pollinator", "Defense"))
})

# Figure 5 in paper. Selection imposed by each of 3 agents: defense, herbivores, pollinators
scaleFUN <- function(x) sprintf("%.2f", x)
plotAgent.Med.Sel <- ggplot(AgentMedSel, aes(x = Agent, y = Gradient, group = Trait))+
  geom_point(size = 3, aes(shape = Trait), position = position_dodge(width = 0.3), alpha = 0.4)+
  xlab("Agent")+ylab("Strength of agent-mediated selection") +
  geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - (1.96*(sd(x)/sqrt(length(x)))),
               fun.ymax = function(x) mean(x) + (1.96*(sd(x)/sqrt(length(x)))),
               geom = "errorbar", width = 0.15, color = "black", aes(group = Agent)) +
  stat_summary(fun.y = mean, geom = "point", size = 5, color = "black", aes(group = Agent)) +
  coord_cartesian(ylim = c(-0.4, 0.55)) + scale_y_continuous(breaks = seq(from = -0.4, to = 0.55, by = 0.1), labels = scaleFUN) +
  ng1 + theme(legend.title=element_blank())
plotAgent.Med.Sel


plotAgent.Med.Sel_ABS <- ggplot(AgentMedSel, aes(x = Agent, y = abs(Gradient), group = Trait))+
  geom_point(size = 3, aes(shape = Trait), position = position_dodge(width = 0.3), alpha = 0.4)+
  xlab("Agent")+ylab("Strength of agent-mediated selection") +
  # geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - (1.96*(sd(x)/sqrt(length(x)))),
               fun.ymax = function(x) mean(x) + (1.96*(sd(x)/sqrt(length(x)))),
               geom = "errorbar", width = 0.15, color = "black", aes(group = Agent)) +
  stat_summary(fun.y = mean, geom = "point", size = 5, color = "black", aes(group = Agent)) +
  coord_cartesian(ylim = c(0, 0.4)) + scale_y_continuous(breaks = seq(from = 0, to = 0.4, by = 0.1), labels = scaleFUN) +
  ng1 + theme(legend.title=element_blank())
plotAgent.Med.Sel_ABS

# Means of absolute value selection gradients
means_AbsVal_Sel <- AgentMedSel %>%
  group_by(Agent) %>%
  summarize(mean = mean(abs(Gradient)))

ggsave("HPS_figures/Figure.5_Sel.x.Agent.pdf", plot = plotAgent.Med.Sel, width = 8, height = 8, unit = "in", dpi = 600)
ggsave("HPS_figures/Figure.5.Inset_ABS_Sel.x.Agent.pdf", plot = plotAgent.Med.Sel_ABS, width = 8, height = 8, unit = "in", dpi = 600)

