if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("GRmetrics")

trace(utils:::unpackPkgZip, edit=TRUE)
install_github("chrisamiller/fishplot")
install.packages('vegan')

library(vegan)
library(ggplot2)
library(devtools)
library(fishplot)
library(agricolae)
library(gridExtra)
library(ggbiplot)
library(GRmetrics)

##Kmeans cluster analysis_Figure S2########################################################################################################
##########################################################################################################
###################################################################################
#########################################################################################
traits.treat<-read.table("Supplementary dataset1_sheet1.txt", header=T, sep="\t", dec=",")
  traits.nomis<- na.omit(traits.treat)
  traits.scale<-scale(traits.nomis[6:19])

  traits14.k<-cascadeKM(traits.scale, inf.gr=2, sup.gr=9, iter=200, criterion="ssi")
  plot(traits14.k, sortg=T)
  #The best partition is indicated by the highest SSI value
  #the 5 groups cluster result was chosed based on the SSI value
  traits14.k$results
  ngroup<-traits14.k$partition

##########Figure 1 ##############################################################################################################################################
########################################################################################################################################################
################################################################################################################################################
#figure 1a
cylce1growth<-read.table("growth reduction after first cycle.txt", header=T, sep="\t", dec=",")
  cylce1growth.lm<-lm(Aboveground.dry.biomass~Treatment, data=cylce1growth)
  summary(cylce1growth.lm)
  boxplot(Aboveground.dry.biomass~Treatment, data=cylce1growth, ylab="Aboveground dry biomass(g)",xlab="", col= 'darkgrey')
  points(Aboveground.dry.biomass~Treatment, data=cylce1growth,col='black',pch=16)
  text(2,0.018,'***')

#figure1 b-f
plant<-read.table("Supplementary dataset2_sheet2.txt", header=T, sep="\t", dec=",")
  plant.nomis<- na.omit(plant)
  names(plant)
  plant.mean=aggregate(.~Sample.ID+Plant.Cycle+Replicate.Line+Phenotypes+Mutations,data=plant.nomis,mean)
  #mean value for each 30 isolates

  mock=plant.mean[which(plant.mean$Phenotypes=="Mock" ),]#extract mock treatment
  plant1=plant.mean[which(!plant.mean$Phenotypes=="Mock" ),]#remove the mock treatment for later plot making

  #plant1$Phenotypes, reorder the different phenotypes 
  plant1$Phenotypes=factor(plant1$Phenotypes)
  plant1$Phenotypes=factor(plant1$Phenotypes, levels = levels(plant1$Phenotypes)[c(1,2,6,5,3,4)]) #choose your order
  names(plant.mean)
  names(plant1)
  col1=c("dimgrey","dimgrey","grey","red","limegreen","darkgreen")#define colors for all phenotypes


####figure 1b
  boxplot(Shoot.biomass.g.~Phenotypes, data=plant1,las = 3, ylab="Shoot biomass (g)", col=col1)
  abline(h=mean(mock$Shoot.biomass.g.), col="black", lwd=2, lty=2)
  RLlm=lm(Shoot.biomass.g.~relevel(Phenotypes,"Mock"), data=plant.mean)
  summary(RLlm)
#extract the significant differences between control plants(Mock) and bacteria-treated plants; ANOVA result

####figure 1c
  boxplot(Root.biomass.g.~Phenotypes, data=plant1,las = 3, ylab="Root biomass (g)", col=col1)
  abline(h=mean(mock$Root_biomass.g.), col="black", lwd=2, lty=2)
  RLlm=lm(Root_biomass.g.~relevel(Phenotypes,"Mock"), data=plant.mean)
  summary(RLlm)

####figure 1d
  boxplot(Root.length.cm.~Phenotypes, data=plant1,las = 3, ylab="Root length (cm)", col=col1)
  abline(h=mean(mock$Root_length.cm.), col="black", lwd=2, lty=2)
  RLlm=lm(Root.length.cm.~relevel(Phenotypes,"Mock"), data=plant.mean)
  summary(RLlm)

####figure 1e
  boxplot(Lateral_root.count.~Phenotypes, data=plant1,las = 3, ylab="Number of lateral roots", col=col1)
  abline(h=mean(mock$Lateral_root.count.), col="black", lwd=2, lty=2)
  RLlm=lm(Lateral_root.count.~relevel(Phenotypes,"Mock"), data=plant.mean)
  summary(RLlm)

####figure 1f
  boxplot(Green.leaf.biomass.pixels.~Phenotypes, data=plant1,las = 3, ylab="Green leaf biomass (pixels)", col=col1)
  abline(h=mean(mock$Green.leaf.biomass.pixels.), col="black", lwd=2, lty=2)
  RLlm=lm(Green.leaf.biomass.pixels.~relevel(Phenotypes,"Mock"), data=plant.mean)
  summary(RLlm)


##########Figure 2#############################################################################
################################################################################
################################################################

####figure2a_Fishplot 
  #base on 'Supplementary dataset1_sheet3': the change of each phenotype proportion
  #sum of five replicates
  parents = c(0,0,0,0,0)
  timepoints=c(0,2,4,6) #plant cycle
  frac.table = matrix(
   c(00,00, 00,00,100,
     00,03, 30,00,67,
     05,13, 21,05,56,
     09,24,00, 41,26),
  ncol=length(timepoints))
  #the value within the matrix are the percentage(rounded up) out of 16*5=80 isolates at each time point
  #create a fish object
  fish = createFishObject(frac.table,parents,timepoints=timepoints)
  #calculate the layout of the drawing
  fish = layoutClones(fish,separate.independent.clones=F)
  fish = setCol(fish,c("red","darkgreen","grey","limegreen","dimgrey"))
  vlines=c(0,2,4,6)
  #draw the plot, using the splining method 
  #and providing both timepoints to label and a plot title
  p0=fishPlot(fish,shape="spline",cex.title=1, vlines=vlines, vlab=vlines, title.btm="Mean", cex.vlab=1,bg.type = "solid", bg.col ="lightyellow")
  ##Draw fishplot for each replicate line 
 
  #e.g. Line1
  parents = c(0,0,0)
  timepoints=c(0,2,4,6)
  frac.table = matrix(
   c(00,00,100,
     00,06,94,
     06,31,63,
     75,00,25),
   ncol=length(timepoints))

  fish = createFishObject(frac.table,parents,timepoints=timepoints)
  fish = layoutClones(fish,separate.independent.clones=F)
  fish = setCol(fish,c( "green","grey", "dimgrey"))
  vlines=c(0,2,4,6)
  p1=fishPlot(fish,shape="spline",cex.title=1, vlines=vlines, vlab=vlines, cex.vlab=1,bg.type = "solid", bg.col ="lightyellow",title.btm="Line 1")

####(figure 2b)
  plant.traits<-scale(plant.mean[,8:11])
  plant.pca<- prcomp(plant.traits)
  groups<-as.factor(plant.mean[,4])
  ggbiplot(plant.pca, obs.scale = 1, var.scale = 1,
         groups = groups, ellipse = TRUE) +
   scale_color_manual(name="Phenotypes", values=c("dimgrey","dimgrey","blue","limegreen","darkgreen","red","grey")) +
   theme_bw()+theme( legend.direction = 'horizontal', legend.position = 'top',panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  plant.performance.index<-(plant.pca $x[,1])*-1

####figure 2c
pminteraction<-read.table("Supplementary dataset2_sheet1.txt", header=T, sep="\t", dec=",")
  names(pminteraction)
  fit=lm(Plant.performance..PC1..1.~Bacterial.population.size..cells.per.plant., data=pminteraction)
  summary(fit)

  col=c("dimgrey","dimgrey","limegreen","darkgreen","red","grey")[pminteraction$Phenotype]
  ggplot(pminteraction, aes(x=Bacterial.population.size..cells.per.plant., y=Plant.performance..PC1..1.)) + 
    geom_point(size=2, color=col)+
    geom_smooth(method='lm',col='black')+
    xlab("Bacterial population (cells per plant)") + ylab("Plant performance(PC1)")+
    theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##########Figure 3#####################################################################################################################
#####################################################################################################################################
####_potential mechanisms of mutualism establishment###########################################################################################################################
 par(mfrow=c(2,2))

####figure 3a
  #firstly PC analysis baterial growth in 14 carbons. 
  #Then extract PC1 value (normalised) as the index of bacterial growth on Arabidopsis root secreted carbon sourses
  carbon<-read.table("Supplementary dataset1_sheet2.txt", header=T, sep="\t", dec=",")
  carbon14<-scale(carbon[,5:18])
  carbon.pca<- prcomp(carbon14)
  groups<-as.factor(carbon[,4])
  ggbiplot(carbon.pca, obs.scale = 1, var.scale = 1,
         groups = groups, ellipse = TRUE) +
   scale_color_manual(name="Phenotypes", values=c("dimgrey","dimgrey","limegreen","darkgreen","red","grey")) +
   theme_bw()+theme( legend.direction = 'horizontal', legend.position = 'top',panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  carbon.index<-(carbon.pca $x[,1])*-1
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  carbon.usage<-range01(carbon.index) #normalise the data from 0 to 1
  
  #write 'carbon.index'as 'carbon pc1' to 'Supplementary dataset1_sheet4'
  #to change to order of different boxws
  carbon$Phenotype=factor(carbon$Phenotype,levels = levels(carbon$Phenotype)[c(1,2,6,5,3,4)])
  boxplot(carbon.usage~Phenotype, data=carbon,las =2,ylab="Growth on Arabidopsis root secreted carbons", col=(c("dimgrey","dimgrey","grey","red","limegreen","darkgreen")))
  carbon.usage.lm<-lm(carbon.usage~Phenotype, data=carbon)
  summary(carbon.usage.lm)
  comparison <- HSD.test(carbon.usage.lm,"Phenotype", group=TRUE,
                       main="Growth on Arabidopsis root secreted carbons")
 
 
####figure 3b
  induced.expression.fold.lm<-lm(induced.expression.fold~Phenotype, data=pminteraction)
  summary(induced.expression.fold.lm)
  comparison <- HSD.test(induced.expression.fold.lm,"Phenotype", group=TRUE,
                         main='Fold induction of MYB72')
  comparison 
  col1=c("dimgrey","dimgrey","grey","red","limegreen","darkgreen")#define colors for all phenotypes
  pminteraction$Phenotype=factor(pminteraction$Phenotype, levels = levels(pminteraction$Phenotype)[c(1,2,6,5,3,4)]) #choose your order
  boxplot(induced.expression.fold~Phenotype, data=pminteraction,las = 3, ylab="Fold change of MYB72 expression",xlab="", col=col1)


####figure 3c
  #data(inputCaseC)
  scopletin<-read.table("Supplementary dataset2_sheet3.txt", header=T, sep="\t", dec=",")
  #Subset the data with starting point and after 72h
  scopletin=subset(scopletin,time=='0' | time=='72' )
  names(scopletin)
  drc_output = GRfit(scopletin, groupingVariables =c('Sample.ID','Phenotype') ,case = "C")
  GRbox(drc_output, metric ='Emax', groupVariable = 'Phenotype', 
        pointColor = 'Sample.ID')
  GRbox(drc_output, metric ='GRmax', groupVariable = 'Phenotype', 
        pointColor = 'Sample.ID')
  GRgetMetrics(drc_output)$Emax
  GRgetMetrics(drc_output)$GRmax
  #extract the metrics, GRmax and Emax, then write it to file 'Supplementary dataset2_sheet1'
  #ANOVA and HSD.test
  #data: Supplementary dataset2_sheet1
  Emax.lm<-lm(Emax~Phenotype, data=pminteraction)
  summary(Emax.lm)
  comparison <- HSD.test(Emax.lm,"Phenotype", group=TRUE,
                         main='Scopletin tolerance')
  comparison 
  boxplot(Emax~Phenotype, data=pminteraction,las = 3, ylab="Scopoletin tolerance (%)",xlab="", col=col1)

####figure 3d
  plot(Emax~induced.expression.fold,data=pminteraction, col=col1[Phenotype],pch=16,ylab="Scopoletin tolerance (%)",xlab="Induced MYB72 expression")
  fit=lm(Emax~induced.expression.fold, data=pminteraction)
  abline(lm(Emax~induced.expression.fold, data=pminteraction),lwd=2)
  summary(fit)
 


########Figure 4#######################################################################################################################
#####################################################################################################################################
####genotype-phenotype###########################################################################################################################
  #data Supplementary dataset2_sheet2: 30 isolates, 3 replicates for each
    
  plant.per<-scale(plant[,8:11])
  plant.pca<- prcomp(plant.per)
  plant.1<-(plant.pca $x[,1])
  # extract PC1 value for every three replicates of each 30 islolate, allowing for later ANOVA analysis
  ## need to scale base on mock treatment value
  plant.performance=cbind(plant, plant.1)
  names(plant.performance)
  plant.mean=aggregate(.~Sample.ID+Plant.Cycle+Replicate.Line+Phenotypes+Mutations,data=plant.performance,mean)
  ggplot(plant.mean, aes(x=Replicate.Line, y=plant.1, colour=Phenotypes))+  
    scale_color_manual(values=c("dimgrey","dimgrey","purple","limegreen","darkgreen","red","grey"))+
    geom_line(col="black") +
    geom_point(size=4)+
    theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  ##statistic analysis for all lines
  #e.g line1
  line1=plant.performance[which(plant.performance$Replicate.Line=='line 1'| plant.performance$Replicate.Line=="ori"),]
  plant1.lm=lm(plant.1~Mutations, data=line1)
  comparison <- HSD.test(plant1.lm,"Mutations", group=TRUE)
  comparison 
  

########Figure 5###########################################################################################################
#####################################################################################################################
###################################################################################################################
##### relative fitness figure for gac mutations
fitness <-read.table("relative fitness gac.txt", header=T, sep="\t", dec=",")
  ggplot(fitness, aes(x=Context, y=r)) + facet_wrap(~mutation, scales = "fix")+geom_hline(yintercept = 1,linetype="dashed")+
    geom_boxplot(fill="grey")+ylab("Relative fitness(r)")+ scale_x_discrete(limits=c("Plant", "KB","LB", "TSB"))+
    theme_bw()+theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  #ANOVA 
  #example
  f=lm(r~Context*mutation, data=fitness)
  summary(f)
  
  f=lm(r~mutation, data=subset(fitness, fitness$Context=="Plant"))
  summary(f)
  comparison <- HSD.test(f,"mutation", group=TRUE)
  comparison 
  
  f=lm(r~Context, data=subset(fitness, fitness$mutation=="gacA(D54Y)"))
  summary(f)
  comparison <- HSD.test(f,"Context", group=TRUE)
  comparison 
  
  
  
#########Figure S3-4######################################################################################################################
#####################################################################################################################################
###############################################################################################################################
  traits.nomis$Phenotype=factor(traits.nomis$Phenotype, levels=levels(traits.nomis$Phenotype)[c(1,2,6,5,3,4)])
  names(traits.nomis)

####figure s3
  #a
  boxplot(Exprotese..OD440.~Phenotype, data=traits.nomis,las = 3,  ylab="Proteolytic activity (OD440) ", col=(c("dimgrey","dimgrey","grey","red","green","darkgreen")))
   f=lm(Exprotese..OD440.~ Phenotype, data=traits.nomis)
   summary(f)
   comparison <- HSD.test(f,"Phenotype", group=TRUE, main='proteolytic activity')
   comparison

   #b
  boxplot(Verticillium.Inhibition....~Phenotype, data=traits.nomis,las = 3,  ylab="Proteolytic activity (OD440) ", col=(c("dimgrey","dimgrey","grey","red","green","darkgreen")))
   f=lm(Verticillium.Inhibition....~ Phenotype, data=traits.nomis)
   summary(f)
   comparison <- HSD.test(f,"Phenotype", group=TRUE,main='antifungal activity')
   comparison

####figure s4
   #a
   boxplot(Biofilm..OD590.~Phenotype, data=traits.nomis,las = 3,  ylab="Biofilm formation (OD590) ", col=(c("dimgrey","dimgrey","grey","red","green","darkgreen")))
   f=lm(Biofilm..OD590.~ Phenotype, data=traits.nomis)
   summary(f)
   comparison <- HSD.test(f,"Phenotype", group=TRUE,main='biofilm formation')
   comparison
   
  #b 
   #firstly PC analysis of ability to grow in the presence of sub lethal doses of the antibiotics streptomycin, tetracycline and penicillin. 
   #Then extract PC1 value (normalised) 
   biostress<-scale(traits.nomis[,15:17])
   bio.pca<- prcomp(biostress)
   groups<-as.factor(traits.nomis[,5])
   ggbiplot(bio.pca, obs.scale = 1, var.scale = 1,
            groups = groups, ellipse = TRUE) +
     scale_color_manual(name="Phenotypes", values=c("dimgrey","dimgrey","grey","red","limegreen","darkgreen")) +
     theme_bw()+theme( legend.direction = 'horizontal', legend.position = 'top',panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
   
   bio.index<-(bio.pca $x[,1])
   range01 <- function(x){(x-min(x))/(max(x)-min(x))}
   bio.index.n=range01(bio.index)
   #write 'carbon.index'as 'carbon pc1' to 'Supplementary dataset1_sheet4'
   boxplot(bio.index.n~Phenotype, data=traits.nomis,las =2,ylab="Biotic stress resistance index", col=(c("dimgrey","dimgrey","grey","red","limegreen","darkgreen")))
   bio.index.n.lm<-lm(bio.index.n~Phenotype, data=traits.nomis)
   summary(bio.index.n)
   comparison <- HSD.test(bio.index.n.lm,"Phenotype", group=TRUE,
                          main="Biotic stress resistance index")
  #c 
   #firstly PC analysis of combined ability of each isolate to grow under oxidative stress, water potential stress and salt stres
   #Then extract PC1 value (normalised) 
   abiostress<-scale(traits.nomis[,c(14,18:19)])
   abio.pca<- prcomp(abiostress)
   groups<-as.factor(traits.nomis[,5])
   ggbiplot(abio.pca, obs.scale = 1, var.scale = 1,
            groups = groups, ellipse = TRUE) +
     scale_color_manual(name="Phenotypes", values=c("dimgrey","dimgrey","grey","red","limegreen","darkgreen")) +
     theme_bw()+theme( legend.direction = 'horizontal', legend.position = 'top',panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
   
   abio.index<-(abio.pca $x[,1])*-1
   range01 <- function(x){(x-min(x))/(max(x)-min(x))}
   abio.index.n=range01(abio.index)
   boxplot(abio.index.n~Phenotype, data=traits.nomis,las =2,ylab="Abiotic stress resistance index", col=(c("dimgrey","dimgrey","grey","red","limegreen","darkgreen")))
   abio.index.n.lm<-lm(abio.index.n~Phenotype, data=traits.nomis)
   summary(abio.index.n)
   comparison <- HSD.test(abio.index.n.lm,"Phenotype", group=TRUE,
                          main="Abiotic stress resistance index")
   
    

###############################################################################################################################
#####################################################################################################################################
####Figure S5_###########################################################################################################################
#pminteraction<-read.table("Supplementary dataset2_sheet1.txt", header=T, sep="\t", dec=",")
    names(pminteraction)
    pminteraction<- na.omit(pminteraction)
    fit=lm(Plant.performance.promotion....~End.proportion.of.phenotype...., data=pminteraction)
    summary(fit)
    ggplot(pminteraction, aes(x=End.proportion.of.phenotype...., y=Plant.performance.promotion...., colour=Phenotype, shape=Replicate.line, fill=Phenotype)) + 
      scale_color_manual(values=c("dimgrey","green","darkgreen","red","grey"))+
      scale_shape_manual(values=c(21, 22, 23,24,25))+
      scale_fill_manual(values=c("dimgrey","green","darkgreen","red","grey"))+
      geom_point(size=2)+
      geom_abline(intercept = 31.845 , slope = 1.636)+
      xlab("Frequency at end time point") + ylab("Effects on plant performance (%)")+
      theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

###############################################################################################################################
#####################################################################################################################################
####Figure Sx_###########################################################################################################################

