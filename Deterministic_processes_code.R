#This is the full code for figures and tables in PAPER NAME
library(phyloseq)
library(ggplot2)


#set file names
com.file<-"EVO_comm_all.tabular"
tree.file<-"EVO_tree_all.nwk"
treat.file<-"EVO_treatment_all_all"
clas.file<-"EVO_clasf_all.tabular"
geo.file<-"EVO_geochem_all.csv"
comm=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))

#community, not transposed
comm2=(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))

#tree
tree=read.tree(file = tree.file)

#taxonomy table
clas=read.table(clas.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)

#treatment table
treat=read.table(treat.file, header = TRUE, sep = "\t", row.names = 1,
                 as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                 check.names = FALSE)

#geochemistries; from file or collaborator's
EVO_all_geochem_fromAP <- 
  #as.data.frame(read.csv("https://www.dropbox.com/s/jtgl97o6w78qqtr/col-reduced_trim_0_EMR_0.2.csv?dl=1", 
  #                     sep = ",", row.names = 1, stringsAsFactors = FALSE ))
  as.data.frame(read.csv(file=geo.file,sep = ",", row.names = 1, stringsAsFactors = FALSE))
#reduce it to be useful to add to treatment table (once 2017_02 phyloseq object is generated)
EVO_all_geochem_fromAP[,-c(1:5, 18), drop=FALSE]->EVO_all_geochem


#need to convert to matrices so phyloseq isn't angry
as.matrix(comm2)->commmatrix
as.matrix(clas)->clasmatrix

#construct the object
ps.EVO_all<-phyloseq(sample_data(treat), 
                     otu_table(commmatrix, taxa_are_rows = TRUE), 
                     tax_table(clasmatrix), 
                     phy_tree(tree))

#check sample sums - remove sums below cutoff (10K?)
sample_sums(ps.EVO_all)
min(sample_sums(ps.EVO_all))

#remove ghost ASVs
prune_taxa(taxa_sums(ps.EVO_all)>0, ps.EVO_all)->ps.EVO_all

#make a function as this will be needed later
ghostbust<-function(physeq){
  prune_taxa(taxa_sums(physeq)>0, physeq)->physeq
}

#rarefy entire dataset
#alternatively, rarefy individual datasets
set.seed(0451)
rarefy_even_depth(ps.EVO_all, 
                  sample.size = min(sample_sums(ps.EVO_all_pruned)), 
                  rngseed =0451, trimOTUs = TRUE)->ps.EVO_all_rar
#this rarefies these data to 10449 reads each

ghostbust(ps.EVO_all_rar)->ps.EVO_all_rar
any(taxa_sums(ps.EVO_all_rar)==0) #just checking



#Because of the random nature of rarefying, write to file
write.table(tax_table(ps.EVO_all_rar),file='EVO_clasf_all_rar.tabular', quote=F, sep='\t')
write.table(otu_table(ps.EVO_all_rar),file='EVO_comm_all_rar.tabular', quote=F, sep='\t')
castor::write_tree(phy_tree(ps.EVO_all_rar),file = "EVO_tree_all_rar.nwk")

#saving to reproduce
save(ps.EVO_all, file = "ps.EVO_all.rda")
save(ps.EVO_all_rar, file = "ps.EVO_all_rar.rda")

subset_samples(ps.EVO_all_rar, FilterSize==0.2)->ps.EVO_all_02_rar

#add the geochems to the 0.2 filter samples
EVO_all_geochem[sample_names(ps.EVO_all_02_rar),,drop=FALSE]->EVO_all_geochem_matched

#merge the sample datas
merge(sample_data(ps.EVO_all_02_rar), EVO_all_geochem_matched, by="row.names", all.x=T)->treat.all
treat.all$Row.names->rownames(treat.all)
treat.all[,-1,drop=FALSE]->treat.all
treat.all->sample_data(ps.EVO_all_02_rar)

#format variables, some to factors, some to timepoints, etc.
as.factor(sample_data(ps.EVO_all_02_rar)$DaysAfter)->sample_data(ps.EVO_all_02_rar)$DaysAfterfactor
as.factor(sample_data(ps.EVO_all_02_rar)$WellID)->sample_data(ps.EVO_all_02_rar)$WellID
as.factor(sample_data(ps.EVO_all_02_rar)$well_type)->sample_data(ps.EVO_all_02_rar)$well_type

treat.all<-data.frame(sample_data(ps.EVO_all_02_rar))
treat.all$timepoint<-recode_factor(treat.all$DaysAfterfactor,
                                   "-6"='1','1'='2','8'='3','15'='4','22'='5',
                                   '50'='6','78'='7','106'='8','134'='9',"-28"='1',
                                   '4'='2','17'='3','31'='4','80'='5','140'='6','269'='7',"372"="10")
treat.all$phase<-recode_factor(treat.all$DaysAfterfactor,
                               "-6"='1','1'='1','8'='2','15'='2','22'='2',
                               '50'='2','78'='3','106'='3','134'='3',"-28"='1',
                               '4'='2','17'='2','31'='2','80'='2','140'='3','269'='3','372'='3')
treat.all$phase->EVOphases
levels(EVOphases)<-c("1",'2','3','Control')

EVOphases[which(treat.all$well_type=="Control")]<-"Control"
EVOphases->treat.all$phase

treat.all->sample_data(ps.EVO_all_02_rar)
ghostbust(ps.EVO_all_02_rar)->ps.EVO_all_02_rar

any(taxa_sums(ps.EVO_all_02_rar)==0)

subset_samples(ps.EVO_all_02_rar, EVOYear==2009)->ps.EVO_2009_rar
subset_samples(ps.EVO_all_02_rar, EVOYear==2017)->ps.EVO_2017_02_rar


#remove ghosts using prune_taxa
prune_taxa(taxa_sums(ps.EVO_2009_rar)>0, ps.EVO_2009_rar)->ps.EVO_2009_rar
prune_taxa(taxa_sums(ps.EVO_2017_02_rar)>0, ps.EVO_2017_02_rar)->ps.EVO_2017_02_rar

#confirm
length(which(taxa_sums(ps.EVO_2009_rar)==0))
length(which(taxa_sums(ps.EVO_2017_02_rar)==0))

#save them
save(ps.EVO_2009_rar, file = "ps.EVO_2009_rar.rda")
save(ps.EVO_2017_02_rar, file = "ps.EVO_2017_02_rar.rda")


###############################################################################
###############################################################################
######Figure 1 & Table 1: Geochemistries
###############################################################################
###############################################################################

#preparing for plotting
treat.all->EVO_sampdata
data.frame(sample_data(ps.EVO_all_02_rar))->EVO_sampdata

EVO_sampdata->EVO_sampdata_forgeochems

pivot_longer(data.frame(EVO_sampdata_forgeochems), cols = c(7:18), names_to = "Variable",
             values_to = "Measurement")->EVO_geochem_long

EVO_geochem_plotting<-EVO_geochem_long %>%
  group_by(Variable,DaysAfter,EVOYear,well_type) %>%
  summarise(mean=mean(Measurement),
            se=plotrix::std.error(Measurement))

#selecting acetate and terminal electron acceptors
unique(EVO_geochem_plotting$Variable)[c(1,4,8,10, 12)]->geochemstoplot

#subsetting them; need to remove anything with day 372 as geochemical measurements weren't taken that day
#preparing data to plot
EVO_geochem_plotting_reduced<-subset(EVO_geochem_plotting, DaysAfter!=372)
EVO_geochem_plotting_reduced<-subset(EVO_geochem_plotting_reduced, Variable %in% geochemstoplot)

as.factor(EVO_geochem_plotting_reduced$EVOYear)->EVO_geochem_plotting_reduced$EVOYear
EVO_geochem_plotting_reduced$upper<-EVO_geochem_plotting_reduced$mean+EVO_geochem_plotting_reduced$se
EVO_geochem_plotting_reduced$lower<-EVO_geochem_plotting_reduced$mean-EVO_geochem_plotting_reduced$se


#plotting

p.EVO_U<-ggplot(data=subset(EVO_geochem_plotting_reduced, Variable=="U238_mgL"), aes(x=DaysAfter,y=mean))+
  geom_line(aes(color=EVOYear, linetype=well_type), linewidth=0.5)+
  geom_ribbon(data=subset(EVO_geochem_plotting_reduced, Variable=="U238_mgL" &well_type=='Monitoring'), 
              aes(ymin= lower, ymax= upper, colour=EVOYear, fill=EVOYear),
              alpha=0.2, linetype=0.75)+
  geom_point(aes(color=EVOYear, shape=well_type), size=3)+
  theme_classic(base_size=12)+
  xlab("Days after injection")+ylab("Uranium (mg/L)")+ 
  scale_color_discrete(name="Injection\nyear")+
  scale_fill_discrete(name="Injection\nyear")+
  scale_linetype_manual(values = c('dashed','solid'),name='Well\ntype')+
  scale_shape_manual(values = c(19,17),name='Well\ntype')+
  theme(text = element_text(size = 12), legend.position = 'none')+scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
p.EVO_U

ggsave('U_nolegend.tiff',
       width=90,
       height = 90,
       units='mm',
       dpi=1000
)

p.EVO_Acetate<-ggplot(data=subset(EVO_geochem_plotting_reduced, Variable=="Acetate_mgL"), aes(x=DaysAfter,y=mean))+
  geom_line(aes(color=EVOYear, linetype=well_type), linewidth=0.5)+
  geom_ribbon(data=subset(EVO_geochem_plotting_reduced, Variable=="Acetate_mgL" &well_type=='Monitoring'), 
              aes(ymin= lower, ymax= upper, colour=EVOYear, fill=EVOYear),
              alpha=0.2, linetype=0.75)+
  geom_point(aes(color=EVOYear, shape=well_type), size=3)+
  theme_classic(base_size=12)+
  xlab("Days after injection")+ylab("Acetate (mg/L)")+ 
  scale_color_discrete(name="Injection\nyear")+
  scale_fill_discrete(name="Injection\nyear")+
  scale_linetype_manual(values = c('dashed','solid'),name='Well\ntype')+
  scale_shape_manual(values = c(19,17),name='Well\ntype')+
  theme(text = element_text(size = 12), legend.position = 'none')+scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
p.EVO_Acetate

ggsave('Acetate_nolegend.tiff',
       width=90,
       height = 90,
       units='mm',
       dpi=1000
)

p.EVO_Iron<-ggplot(data=subset(EVO_geochem_plotting_reduced, Variable=="Fe54_mgL"), aes(x=DaysAfter,y=mean))+
  geom_line(aes(color=EVOYear, linetype=well_type), linewidth=0.5)+
  geom_ribbon(data=subset(EVO_geochem_plotting_reduced, Variable=="Fe54_mgL" &well_type=='Monitoring'), 
              aes(ymin= lower, ymax= upper, colour=EVOYear, fill=EVOYear),
              alpha=0.2, linetype=0.75)+
  geom_point(aes(color=EVOYear, shape=well_type), size=3)+
  theme_classic(base_size=12)+
  xlab("Days after injection")+ylab("Iron (mg/L)")+ 
  scale_color_discrete(name="Injection\nyear")+
  scale_fill_discrete(name="Injection\nyear")+
  scale_linetype_manual(values = c('dashed','solid'),name='Well\ntype')+
  scale_shape_manual(values = c(19,17),name='Well\ntype')+
  theme(text = element_text(size = 12),legend.position = 'none' )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

p.EVO_Iron


p.EVO_NO3<-ggplot(data=subset(EVO_geochem_plotting_reduced, Variable=="NO3_mgL"), aes(x=DaysAfter,y=mean))+
  geom_line(aes(color=EVOYear, linetype=well_type), linewidth=0.5)+
  geom_ribbon(data=subset(EVO_geochem_plotting_reduced, Variable=="NO3_mgL" &well_type=='Monitoring'), 
              aes(ymin= lower, ymax= upper, colour=EVOYear, fill=EVOYear),
              alpha=0.2, linetype=0.75)+
  geom_point(aes(color=EVOYear, shape=well_type), size=3)+
  theme_classic(base_size=12)+
  xlab("Days after injection")+ylab("Nitrate (mg/L)")+ 
  scale_color_discrete(name="Injection\nyear")+
  scale_fill_discrete(name="Injection\nyear")+
  scale_linetype_manual(values = c('dashed','solid'),name='Well\ntype')+
  scale_shape_manual(values = c(19,17),name='Well\ntype')+
  theme(text = element_text(size = 12), legend.position = 'none')+scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

p.EVO_NO3

p.EVO_SO4<-ggplot(data=subset(EVO_geochem_plotting_reduced, Variable=="SO4_mgL"), aes(x=DaysAfter,y=mean))+
  geom_line(aes(color=EVOYear, linetype=well_type), linewidth=0.5)+
  geom_ribbon(data=subset(EVO_geochem_plotting_reduced, Variable=="SO4_mgL" &well_type=='Monitoring'), 
              aes(ymin= lower, ymax= upper, colour=EVOYear, fill=EVOYear),
              alpha=0.2, linetype=0.75)+
  geom_point(aes(color=EVOYear, shape=well_type), size=3)+
  theme_classic(base_size=12)+
  xlab("Days after injection")+ylab("Sulfate (mg/L)")+ 
  scale_color_discrete(name="Injection\nyear")+
  scale_fill_discrete(name="Injection\nyear")+
  scale_linetype_manual(values = c('dashed','solid'),name='Well\ntype')+
  scale_shape_manual(values = c(19,17),name='Well\ntype')+
  theme(text = element_text(size = 12), legend.position = 'none')+scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

p.EVO_SO4


#next: spearman's rank correlations between geochemical variables.
#Not necessarily trying to make powerful statistical insights 
#mostly just trying to show directionality of geochems in relation to each other

EVO_sampdata[,c(1:18),drop=FALSE]->cortestdata
#removing timepoint w/ no geochem observations
subset(cortestdata, DaysAfter!=372)->cortestdata

#subsetting by year and running a spearman's rank
subset(cortestdata, EVOYear==2009)->cortestdata09
subset(cortestdata, EVOYear==2017)->cortestdata17
geospearman09<-cor(cortestdata09[,7:18], method="spearman", use="pairwise.complete.obs")
geospearman17<-cor(cortestdata17[,7:18], method="spearman", use="pairwise.complete.obs")

#for corrplotting, putting the 2009 observations in the upper triangle and the 2017 in the lower
geospearman09->geospearman09_forcombn
geospearman09_forcombn[lower.tri(geospearman09_forcombn)]<-geospearman17[lower.tri(geospearman17)]

#making an initial plot
corrplot::corrplot(geospearman09_forcombn, type="full",
                   diag=FALSE, 
                   method="square", 
                   addCoef.col = "black", 
                   tl.pos = "d", 
                   tl.col = 'black', 
                   number.cex=0.7)

#colnames(geospearman09_forcombn)<-c('pH','Specific\nconductivity', 'Mg',
#                                   'Al','K','Ca','Fe','Mn','U','NO\_3','SO4','Acetate')
#rownames(geospearman09_forcombn)<-colnames(geospearman09_forcombn)

#using ggcorrplot for more plotting flexibility
geospearman_corrplot<-ggcorrplot::ggcorrplot(geospearman09_forcombn, method='square',type='full', show.diag = F,
                                             hc.order=T, colors = c('darkblue','white','red'),
                                             legend.title='Spearman\'s\nr')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=0.9),
        axis.title = element_blank(),
        text=element_text(size=12))+
  scale_x_discrete(labels=c('pH',expression('NO'[3]),'U',expression('SO'[4]),
                            expression(kappa), 'Fe','Mn','Acetate','Al','K','Mg','Ca'))+
  scale_y_discrete(labels=c('pH',expression('NO'[3]),'U',expression('SO'[4]),
                            expression(kappa), 'Fe','Mn','Acetate','Al','K','Mg','Ca'))

geospearman_corrplot



#need a function to save legends for figure editing
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

get_legend(p.EVO_U+theme(legend.position = 'top'))->p.EVO_U_legend

ggsave('U_legend.tiff',
       plot=p.EVO_U_legend,
       dpi=1000,
       width=190,
       height=50, units='mm')

Fig_1<-gridExtra::grid.arrange(p.EVO_U+annotate('text',x=max(layer_scales(p.EVO_U)$x$get_limits()),
                                                y=max(layer_scales(p.EVO_U)$y$get_limits()),label='a', size=6),
                               p.EVO_NO3+annotate('text',x=max(layer_scales(p.EVO_NO3)$x$get_limits()),
                                                  y=max(layer_scales(p.EVO_NO3)$y$get_limits()),label='b', size=6),
                               
                               p.EVO_Iron+annotate('text',x=max(layer_scales(p.EVO_Iron)$x$get_limits()),
                                                   y=max(layer_scales(p.EVO_Iron)$y$get_limits()),label='c', size=6),
                               p.EVO_SO4+annotate('text',x=max(layer_scales(p.EVO_SO4)$x$get_limits()),
                                                  y=max(layer_scales(p.EVO_SO4)$y$get_limits()),label='d', size=6),
                               
                               p.EVO_Acetate+annotate('text',x=max(layer_scales(p.EVO_Acetate)$x$get_limits()),
                                                      y=max(layer_scales(p.EVO_Acetate)$y$get_limits()),label='e', size=6),
                               geospearman_corrplot+theme(legend.position = 'left',legend.key.width = unit(3, 'mm'),
                                                          legend.key.height = unit(4, 'mm'),
                                                          legend.title = element_text(size = 6),
                                                          legend.text = element_text(size=6))+
                                 coord_cartesian(clip='off')+
                                 annotate('text', x=12,y=12.5,label='f',size=6),
                               nrow=3, ncol=2)

Fig_1

ggsave('Fig1.tiff',
       plot = Fig_1,
       width=190,
       height = 240,
       units='mm',
       dpi=1000
)
export::graph2ppt()


#mantel test to determine how similar these matrices really are
vegan::mantel(geospearman09, geospearman17, method = "spearman", permutations = 999)


###############################################################################
###############################################################################
###### Figure 2, Figures S1 and S2: Abundance and alpha diversity
###############################################################################
###############################################################################
###### Abundance


#importing abundance data from Gihring et al, 2011 as well as the 2017 injection
library(readxl)
AODC_EVO <- read_excel("~/R_files/EVO_09vs17_inputs/AODC_EVO.xlsx", 
                       col_types = c("text", "skip", "numeric", 
                                     "numeric", "text"))
Gihring_Abundance <- read.csv("~/R_files/EVO_09vs17_inputs/Gihring Abundance.csv")
Gihring_Abundance[-7, ,drop=FALSE]->Gihring_Abundance_reduced
#changing to LFC
Gihring_Abundance_reduced$logfoldchange<-log2(Gihring_Abundance_reduced$Normtot)
Gihring_Abundance_reduced[,c(1:2,15),drop=FALSE]->abundancelfc_2009
abundancelfc_2009$WellID<-factor(abundancelfc_2009$WellID, levels=c("FW215","FW202","MLSG4"))


#further abundances
AODC_EVO <- read_excel("~/R_files/EVO_09vs17_inputs/AODC_EVO.xlsx", 
                       col_types = c("text", "skip", "numeric", 
                                     "numeric", "text"))
AODC_EVO->AODC_EVO_01and02

#replacing TNTC with the max 10^ value in the dataset


AODC_EVO[which((AODC_EVO[,4])=="TNTC"),4]<-"10000000"
AODC_EVO$`AODC_cell/mL`<-as.numeric(AODC_EVO$`AODC_cell/mL`)

#AODC_EVO[15,4]<-10000000
#AODC_EVO[14,4]<-10000000
subset(AODC_EVO, FilterSize==0.2)->AODC_EVO
#subset(AODC_EVO, DaysAfter<=135)->AODC_EVO

#converting to fold-change
AODC_EVO$foldchange<-NA

i=1
for (i in 1:nrow(AODC_EVO)){
  AODC_EVO[i,1]->well
  pull(well)->well
  AODC_EVO[i,5]<-(AODC_EVO[i,4]/AODC_EVO[which(AODC_EVO$WellID %in% well & AODC_EVO$DaysAfter==-6),4])
}
AODC_EVO$logfoldchange<-log2(AODC_EVO$foldchange)

#plotting
abundancelfc_2017<-AODC_EVO

#standardize colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
lfc_cols<-c("FW202"=cbPalette[2], 'FW215'=cbPalette[1],'MLSG4'=cbPalette[3],
            'MLSB3'=cbPalette[4],'FW216'=cbPalette[5], 'GP01'=cbPalette[6],'GP03'=cbPalette[7])

p_abundancelfc_2017_bars<-ggplot(data=subset(abundancelfc_2017,DaysAfter!=-6), 
                                 aes(x=as.factor(DaysAfter), y=logfoldchange))+
  geom_bar(stat="identity",position="dodge", aes(fill=WellID),width = 0.5)+theme_classic(base_size=15)+
  guides(fill=guide_legend(title="Well ID"))+
  ylab("Log2-fold change in\nAODC from preinjection") + xlab("Days after injection")+
  geom_hline(yintercept=0)+ggtitle("b")+  
  theme(text=element_text(size=14,  family="sans"),plot.title=element_text(hjust=1.25))+
  scale_fill_manual(values=lfc_cols)


p_abundancelfc_2009_bars<-ggplot(data=subset(abundancelfc_2009, DaysAfter !=-28), 
                                 aes(x=as.factor(DaysAfter), y=logfoldchange))+
  geom_bar(stat="identity",position="dodge", aes(fill=WellID),width=0.5)+theme_classic(base_size=15)+
  ylab("Log2-fold change in\n16S rRNA gene copies/L") + xlab("Days after injection")+
  geom_hline(yintercept=0)+ggtitle("a") + 
  theme(text=element_text(size=14,  family="sans"),plot.title=element_text(hjust=1))+
  guides(fill=guide_legend(title="Well ID"))+ggtitle("a") + theme(plot.title=element_text(hjust=1.25))+
  scale_fill_manual(values=lfc_cols)


Fig_S1<-gridExtra::grid.arrange(p_abundancelfc_2009_bars, p_abundancelfc_2017_bars)
ggsave('Figure_S1.tiff',
       plot = Fig_S1,
       width=190,
       height = 190,
       units='mm',
       dpi=1000
)


###############################################################################
###### alpha diversity

#this is our sample data for every sample
EVO_sampdata 

#set a new column with the negative values converted to 0 to avoid some collisions down the line
EVO_sampdata$DaysAfterfactor_0<-as.factor(EVO_sampdata$DaysAfter)
EVO_sampdata$DaysAfterfactor_0<-dplyr::recode_factor(EVO_sampdata$DaysAfterfactor,
                                                     "-28"="0", "-6"="0")
#libraries for alpha diversity metrics
library(ieggr)
library(picante)


alpha.g(t(otu_table(ps.EVO_all_02_rar)), 
        td.method = c("richness", "shannon", "simpson"))->EVO_all_rar_alpha

picante::pd(samp =t(otu_table(ps.EVO_all_02_rar)) , 
            tree = phy_tree(ps.EVO_all_02_rar),
            include.root = FALSE)->EVO_all_rar_palpha

merge(EVO_all_rar_palpha, EVO_all_rar_alpha, by=0)->EVO_all_rar_alpha
EVO_all_rar_alpha[,1]->rownames(EVO_all_rar_alpha)
EVO_all_rar_alpha[,-1]->EVO_all_rar_alpha
EVO_all_rar_alpha[,-2]->EVO_all_rar_alpha

#trying some evenness metrics:
library(microbiome)
#tabula package masks evenness 
microbiome::evenness(ps.EVO_all_02_rar)->EVO_evenness_all
microbiome::dominance(ps.EVO_all_02_rar, index = "all", relative=TRUE)->EVO_dominance_all
merge(EVO_all_rar_alpha, EVO_evenness_all, by=0)->EVO_all_rar_alpha
rownames(EVO_all_rar_alpha)<-EVO_all_rar_alpha$Row.names
EVO_all_rar_alpha[,-1,drop=FALSE]->EVO_all_rar_alpha
merge(EVO_all_rar_alpha, EVO_dominance_all, by=0)->EVO_all_rar_alpha
rownames(EVO_all_rar_alpha)<-EVO_all_rar_alpha$Row.names
EVO_all_rar_alpha[,-1,drop=FALSE]->EVO_all_rar_alpha
merge(EVO_sampdata, EVO_all_rar_alpha, by=0)->EVO_sampdata_all
rownames(EVO_sampdata_all)<-EVO_sampdata_all$Row.names
EVO_sampdata_all[,-1]->EVO_sampdata_all
#need to update EVO_sampdata_all with the pipeline file
merge(EVO_sampdata_all, EVO_allenh, by=0)->EVO_sampdata_all
rownames(EVO_sampdata_all)<-EVO_sampdata_all$Row.names
EVO_sampdata_all[,-1]->EVO_sampdata_all

#adding it back to the sample data object
EVO_sampdata_all->sample_data(ps.EVO_all_02_rar)
###############################################################################
#first, shannon
library(multcompView)
#need this function:

#functions adapted from amplicon::alpha_boxplot

generate_label_df = function(TUKEY, variable) {
  Tukey.levels = TUKEY[[variable]][, 4]
  Tukey.labels = data.frame(multcompView::multcompLetters(Tukey.levels)["Letters"])
  Tukey.labels$group = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$group), 
  ]
  return(Tukey.labels)
}

div_comparisons<-function(metadata,alpha_div, groupID, index){
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx, , drop = F]
  alpha_div = alpha_div[rownames(metadata), ]
  sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
  df = cbind(alpha_div[rownames(sampFile), index], sampFile)
  colnames(df) = c(index, "group")
  model = aov(df[[index]] ~ group, data = df)
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  Tukey_HSD_table = as.data.frame(Tukey_HSD$group)
  
  LABELS = generate_label_df(Tukey_HSD, "group")
  df$stat = LABELS[as.character(df$group), ]$Letters
  
  max = max(df[, c(index)])
  min = min(df[, index])
  x = df[, c("group", index)]
  #col=c('~',index)
  y = x %>% group_by(group) %>% summarise(Max = max(!!sym(index)))
  y = as.data.frame(y)
  rownames(y) = y$group
  df$y = y[as.character(df$group), ]$Max + (max - min) * 0.05
  return(df)
}

subset(EVO_sampdata_all, EVOYear=="2017")->EVO_sampdata_17
subset(EVO_sampdata_all, EVOYear=="2009")->EVO_sampdata_09


df09.shannon<-div_comparisons(metadata = subset(EVO_sampdata_09,WellID !="FW215" ),
                              alpha_div=subset(EVO_sampdata_09,WellID !="FW215" ),
                              groupID='DaysAfterfactor_0',
                              index='shannon')
df09.AllenH<-div_comparisons(metadata = subset(EVO_sampdata_09,WellID !="FW215" ),
                             alpha_div=subset(EVO_sampdata_09,WellID !="FW215" ),
                             groupID='DaysAfterfactor_0',
                             index='Allen.H')
df09.faith<-div_comparisons(metadata = subset(EVO_sampdata_09,WellID !="FW215" ),
                            alpha_div=subset(EVO_sampdata_09,WellID !="FW215" ),
                            groupID='DaysAfterfactor_0',
                            index='PD')
df09.richness<-div_comparisons(metadata = subset(EVO_sampdata_09,WellID !="FW215" ),
                               alpha_div=subset(EVO_sampdata_09,WellID !="FW215" ),
                               groupID='DaysAfterfactor_0',
                               index='richness')


p.shannon.09<-ggplot(df09.shannon, aes(x = group, y = shannon)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.5, width = 0.5, fill = "transparent") + 
  
  theme_classic(base_size=15) + geom_text(data = df09.shannon, aes(x = group, 
                                                                   y = y, label = stat)) +
  theme(text = element_text( size = 12))+ geom_point(data=ctrlpoints09,
                                                     shape=16, 
                                                     aes(x=DaysAfterfactor_0, y=shannon), 
                                                     color="red", size=3)+
  xlab("Days after injection") + ylab("Shannon") +
  scale_x_discrete(labels=c("0"="-28"))+theme(text = element_text(size = 20))+ylim(0,8)
p.shannon.09


# 17 shannon

df17.shannon<-div_comparisons(metadata = subset(EVO_sampdata_17,WellID !="FW215" ),
                              alpha_div=subset(EVO_sampdata_17,WellID !="FW215" ),
                              groupID='DaysAfterfactor_0',
                              index='shannon')
df17.AllenH<-div_comparisons(metadata = subset(EVO_sampdata_17,WellID !="FW215" ),
                             alpha_div=subset(EVO_sampdata_17,WellID !="FW215" ),
                             groupID='DaysAfterfactor_0',
                             index='Allen.H')
df17.faith<-div_comparisons(metadata = subset(EVO_sampdata_17,WellID !="FW215" ),
                            alpha_div=subset(EVO_sampdata_17,WellID !="FW215" ),
                            groupID='DaysAfterfactor_0',
                            index='PD')
df17.richness<-div_comparisons(metadata = subset(EVO_sampdata_17,WellID !="FW215" ),
                               alpha_div=subset(EVO_sampdata_17,WellID !="FW215" ),
                               groupID='DaysAfterfactor_0',
                               index='richness')


p.shannon.17<-ggplot(df17.shannon, aes(x = group, y = df17.shannon[["shannon"]])) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.5, width = 0.5, fill = "transparent") + 
  labs(x = "Groups", y = paste("Shannon", "index")) + 
  theme_classic(base_size=15) + geom_text(data = df17.shannon, aes(x = group, 
                                                                   y = y, label = stat)) +
  theme(text = element_text( size = 7))+ 
  scale_x_discrete(labels=c("0"="-6"))+ 
  geom_point(data=ctrlpoints17, shape=16, aes(x=DaysAfterfactor_0, y=shannon), 
             color="red",
             size=3)+
  xlab("Days after injection") + ylab("Shannon")+theme(text = element_text(size = 20))+ylim(0,8)


p.shannon.17

################################################################################
#richness 

p.richness.09<-ggplot(df09.richness, aes(x = group, y = richness)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.5, width = 0.5, fill = "transparent") + 
  labs(x = "Groups", y = paste(index, "index")) + 
  theme_classic(base_size=15) + geom_text(data = df09.richness, aes(x = group, 
                                                                    y = y, label = stat), nudge_y=100) +
  theme(text = element_text(family = "sans", size = 12))+ geom_point(data=ctrlpoints09,
                                                                     shape=16, 
                                                                     aes(x=DaysAfterfactor_0, y=richness), 
                                                                     color="red", size=3)+
  xlab("Days after injection") + ylab("Richness") +
  scale_x_discrete(labels=c("0"="-28"))+theme(text = element_text(size = 20))+ylim(0,2200)
p.richness.09


p.richness.17<-ggplot(df17.richness, aes(x = group, y = richness)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.7, width = 0.5, fill = "transparent") + 
  labs(x = "Groups", y = paste(index, "index")) + 
  theme_classic(base_size=15) + geom_text(data = df17.richness, aes(x = group, 
                                                                    y = y, label = stat)) +
  geom_point(data=ctrlpoints17,
             shape=16, size=3, 
             aes(x=DaysAfterfactor_0, y=richness), 
             color="red")+
  xlab("Days after injection") + ylab("Richness") +scale_x_discrete(labels=c("0"="-6"))+
  theme(text = element_text(size = 20))

p.richness.17



################################################################################

#faith's PD 
#listed as PD in the big sheet


p.faith.09<-ggplot(df09.faith, aes(x = group, y = PD)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.5, width = 0.5, fill = "transparent") + 
  labs(x = "Groups", y = paste("Faith", "index")) + 
  theme_classic(base_size=15) + geom_text(aes(x = group, y = y, label = stat)) +
  geom_point(data=ctrlpoints09,
             shape=16, aes(x=DaysAfterfactor_0, y=PD), color="red",
             size=3)+
  xlab("Days after injection") + ylab("Faith's PD")+ 
  scale_x_discrete(labels=c("0"="-28"))+theme(text = element_text(size = 20))

p.faith.09



#2017


p.faith.17<-ggplot(df17.faith, aes(x = group, y = PD)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.7, width = 0.5, fill = "transparent") + 
  labs(x = "Groups", y = paste(index, "index")) + 
  theme_classic(base_size=15) + geom_text(aes(x = group, y = y, label = stat)) +
  geom_point(data=ctrlpoints17,
             shape=16, 
             aes(x=DaysAfterfactor_0, y=PD), 
             color="red",
             size=3)+
  xlab("Days after injection") + ylab("Faith's PD")+scale_x_discrete(labels=c("0"="-6"))+theme(text = element_text(size = 20))

p.faith.17


#allen's H (phylogenetic shannon)
# 09

p.allen.09<-ggplot(df09.AllenH, aes(x = group, y = Allen.H)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.7, width = 0.5, fill = "transparent") + 
  labs(x = "Groups", y = "Allen's H") + 
  theme_classic(base_size=15) + geom_text(aes(x = group, y = y, label = stat)) +
  theme(text = element_text(family = "sans", size = 12))+ geom_point(data=ctrlpoints09,
                                                                     shape=16,size=3, 
                                                                     aes(x=DaysAfterfactor_0, y=Allen.H), 
                                                                     color="red")+
  xlab("Days after injection") + ylab("Allen's H") +scale_x_discrete(labels=c("0"="-28"))+theme(text = element_text(size = 20))
p.allen.09
export::graph2ppt(file='allen09')


# 17 Allen.H

p.allen.17<-ggplot(df17.AllenH, aes(x = group, y = Allen.H)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.7, width = 0.5, fill = "transparent") + 
  labs(x = "Days after injection", y = "Allen's H") + 
  theme_classic(base_size=15) + geom_text(aes(x = group, y = y, label = stat)) +
  geom_point(data=ctrlpoints17,
             shape=16,size=3, 
             aes(x=DaysAfterfactor_0, y=Allen.H), 
             color="red")+
  xlab("Days after injection") + ylab("Allen's H") +scale_x_discrete(labels=c("0"="-6"))+theme(text = element_text(size = 20))


p.allen.17
export::graph2ppt(file='allen17')

#can also plot evenness metrics, but a little sloppily
evenness<-EVO_sampdata_all[,c("EVOYear","DaysAfterfactor", "pielou")]

df17.pielou<-div_comparisons(metadata = subset(EVO_sampdata_17,WellID !="FW215" ),
                             alpha_div=subset(EVO_sampdata_17,WellID !="FW215" ),
                             groupID='DaysAfterfactor_0',
                             index='pielou')

p.pielou.17<-ggplot(df17.pielou, aes(x = group, y = pielou)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.7, width = 0.5, fill = "transparent") + 
  labs(x = "Days after injection", y = "Allen's H") + 
  theme_classic(base_size=15) + geom_text(aes(x = group, y = y, label = stat)) +
  theme(text = element_text(family = "sans", size = 12))+ geom_point(data=ctrlpoints17,
                                                                     shape=16,size=3, 
                                                                     aes(x=DaysAfterfactor_0, y=pielou), 
                                                                     color="red")+
  xlab("Days after injection") + ylab("Pielou's evenness") +scale_x_discrete(labels=c("0"="-6"))+theme(text = element_text(size = 15))





p.evenness<-ggplot(evenness, aes(x = DaysAfterfactor, y = pielou)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, outlier.size = 0, 
               size = 0.7, width = 0.5, fill = "transparent") + 
  labs(x = "Days after injection", y = paste(index, "index")) + 
  facet_grid(.~EVOYear)+
  theme_classic(base_size=15) +
  theme(text = element_text(family = "sans", size = 12))

#generating figure 2 and figure S2
Fig_2<-gridExtra::grid.arrange(p.shannon.09+ylim(0,9)+
                                 coord_cartesian(clip='off')+
                                 annotate('text',x=8, y=9, label='a',size=6)+
                                 theme(text = element_text(size = 12)),
                               p.shannon.17+ylim(0,9)+
                                 coord_cartesian(clip='off')+
                                 annotate('text',x=10, y=9, label='b',size=6)+
                                 theme(text = element_text(size = 12)),
                               p.faith.09+ylim(0,210)+
                                 coord_cartesian(clip='off')+
                                 annotate('text',x=8, y=210, label='c',size=6)+
                                 theme(text = element_text(size = 12)),
                               p.faith.17+ylim(0,210)+
                                 coord_cartesian(clip='off')+
                                 annotate('text',x=10, y=210, label='d',size=6)+
                                 theme(text = element_text(size = 12)))
ggsave('Figure_2.tiff',
       plot=Fig_2,
       width=190,
       height=190,
       units='mm',
       dpi=1000)

Fig_S2<-gridExtra::grid.arrange(p.richness.09+ylim(0,2800)+
                                  theme(text = element_text(size = 20))+
                                  coord_cartesian(clip='off')+
                                  annotate('text',x=8, y=2800, label='a',size=6),
                                p.richness.17+ylim(0,2800)+
                                  theme(text = element_text(size = 20))+
                                  coord_cartesian(clip='off')+
                                  annotate('text',x=11, y=2800, label='b',size=6),
                                p.allen.09+ylim(0,4)+
                                  theme(text = element_text(size = 20))+
                                  coord_cartesian(clip='off')+
                                  annotate('text',x=8, y=4, label='c',size=6),
                                p.allen.17+ylim(0,4)+
                                  theme(text = element_text(size = 20))+
                                  coord_cartesian(clip='off')+
                                  annotate('text',x=11, y=4, label='d',size=6))
ggsave(filename = 'Figure_S2.tiff',
       plot=Fig_S2,
       device='tiff',
       width=190,
       height=190,
       units='mm',
       dpi=1000)

################################################################################
################################################################################
#######Fig.3 Beta Diversity#############################
################################################################################
################################################################################
library(vegan)
library(ggordiplots)

#need the function from ggordiplots:
scale_arrow <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  
  if (rev[1] < 0) {
    u[1:2] <- u[2:1]
  }
  if (rev[2] < 0) {
    u[3:4] <- u[4:3]
  }
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  invisible(fill * min(u))
}

##############start with 2009
#ordinating with CCA and NMDS
ordinate(ps.EVO_2009_rar, method = "NMDS", distance = "bray")->ord.EVO_09_rar_NMDS_bray
ordinate(ps.EVO_2009_rar, method = "CCA", distance = "bray")->ord.EVO_09_rar_CCA_bray

#using plot_ordination to determine most useful ordination method
plot_ordination(ps.EVO_2009_rar, ordination=ord.EVO_09_rar_CCA_bray)

#running envfit
envfit(ord.EVO_09_rar_CCA_bray, env=sample_data(ps.EVO_2009_rar)[,c(13,15:18)],
       permutations=999, na.rm=TRUE)->envfit.EVO_09_rar_CCA_bray

#preparing for plotting
ord<-ord.EVO_09_rar_CCA_bray
choices<-c(1,2)
scaling<-1
pt.size<-3
angle<-20
len<-0.5
arrow.col<-"black"
unit<-"cm"
df_ord <- vegan::scores(ord, display = "sites", choices = c(1,2), 
                        scaling = 1)
df_ord <- as.data.frame(df_ord)
axis.labels <- ggordiplots::ord_labels(ord)[choices]
#there's some code here depending on "groups", which may be specified in envfit
#a formula may also be specified

#changing to the names used for plotting
colnames(df_ord) <- c("x", "y")
axis.labels <- ggordiplots::ord_labels(ord)[choices]
fit<-envfit.EVO_09_rar_CCA_bray
df_arrows <- as.data.frame(vegan::scores(fit, "vectors"))
mult <- scale_arrow(df_arrows, df_ord[, c("x", "y")])
df_arrows <- mult * df_arrows
df_arrows$var <- rownames(df_arrows)
df_arrows$p.val <- fit$vectors$pvals
colnames(df_arrows) <- c("x", "y", "var", "p.val")
df_arrows <- df_arrows[df_arrows$p.val <= 0.05, ]
df_arrows->df_arrows_09

merge(as.data.frame(sample_data(ps.EVO_2009_rar)), df_ord, by=0)->df_ord_09

p.ordenv.EVO_2009_rar_CCA_bray<-ggplot(data = df_ord_09, aes(x = x, y = y, color = DaysAfterfactor)) + 
  geom_point(size = 2, aes(shape=well_type)) + xlab("CA1 [7.8%]") + ylab("CA2 [7.2%]")+ 
  geom_segment(data = df_arrows_09, aes(x = 0,  xend = x, y = 0, yend = y), 
               arrow = arrow(angle = angle,  length = unit(len, unit)), color = arrow.col) + 
  geom_text(data = df_arrows_09, aes(x = x, y = y, label = var), 
            color = arrow.col, hjust = "outward") 


fig3_a<-p.ordenv.EVO_2009_rar_CCA_bray +theme_classic(base_size=15)+
  scale_color_viridis_d(name="Days after\ninjection",option="D")+
  scale_shape_discrete(name="Well type")

###########################################################################

##########2017
#workflow modified a bit due to lack of env data on day 372


env17<-as.data.frame(sample_data(ps.EVO_2017_02_rar_nod372)[,c(13,15:18), drop=FALSE])

ordinate(ps.EVO_2017_02_rar, method = "NMDS", distance = "bray")->ord.EVO_17_02_rar_NMDS_bray
subset_samples(ps.EVO_2017_02_rar, DaysAfter!=372)->ps.EVO_2017_02_rar_nod372
ordinate(ps.EVO_2017_02_rar_nod372, method = "NMDS", distance = "bray")->ord.EVO_17_02_rar_nod372_NMDS_bray

ordinate(ps.EVO_2017_02_rar, method = "NMDS", distance = "bray")->ord.EVO_17_rar_NMDS_bray
ordinate(subset_samples(ps.EVO_2017_02_rar, 
                        sample_names(ps.EVO_2017_02_rar) %in% rownames(env17)), 
         method = "CCA", distance = "bray")->ord.EVO_17_rar_CCA_bray_forenvfit
ordinate(ps.EVO_2017_02_rar, method = "CCA", distance = "bray")->ord.EVO_17_rar_CCA_bray

plot_ordination(ps.EVO_2017_02_rar, ordination=ord.EVO_17_rar_CCA_bray)


envfit(ord.EVO_17_02_rar_nod372_NMDS_bray, env=env17, na.rm=TRUE,
       permutations=999)->envfit.EVO_17_02_rar_nod372_NMDS_bray

envfit(ord.EVO_17_rar_CCA_bray_forenvfit, env=env17,
       permutations=999, na.rm=TRUE)->envfit.EVO_17_rar_CCA_bray
ord<-ord.EVO_17_rar_CCA_bray
choices<-c(1,2)
scaling<-1
pt.size<-3
angle<-20
len<-0.5
arrow.col<-"black"
unit<-"cm"
df_ord <- vegan::scores(ord, display = "sites", choices = choices, 
                        scaling = scaling)
df_ord <- as.data.frame(df_ord)
axis.labels <- ggordiplots::ord_labels(ord)[choices]
#there's some code here depending on "groups", which may be specified in envfit
#a formula may also be specified

#changing to the names used for plotting
colnames(df_ord) <- c("x", "y")
axis.labels <- ggordiplots::ord_labels(ord)[choices]
fit<-envfit.EVO_17_rar_CCA_bray
df_arrows <- as.data.frame(vegan::scores(fit, "vectors"))
mult <- scale_arrow(df_arrows, df_ord[, c("x", "y")])
df_arrows <- mult * df_arrows
df_arrows$var <- rownames(df_arrows)
df_arrows$p.val <- fit$vectors$pvals
colnames(df_arrows) <- c("x", "y", "var", "p.val")
df_arrows <- df_arrows[df_arrows$p.val <= 0.05, ]
df_arrows->df_arrows_17

#values from plot_ordination function
xlab <- "CA1 [8.9%]"
ylab <- "CA2 [7.4%]"
merge(as.data.frame(sample_data(ps.EVO_2017_02_rar)), df_ord, by=0)->df_ord_17

p.ordenv.EVO_2017_02_rar_CCA_bray<-ggplot(data = df_ord_17, aes(x = x, y = y, color = DaysAfterfactor)) + 
  geom_point(size = 2, aes(shape=well_type)) + xlab(xlab) + ylab(ylab)


p.ordenv.EVO_2017_02_rar_CCA_bray<-p.ordenv.EVO_2017_02_rar_NMDS_bray + 
  geom_segment(data = df_arrows_17, aes(x = 0,  xend = x, y = 0, yend = y), 
               arrow = arrow(angle = angle,  length = unit(len, unit)), color = arrow.col) + 
  geom_text(data = df_arrows_17, aes(x = x, y = y, label = var), 
            color = arrow.col, hjust = "outward")  

fig3_b<-p.ordenv.EVO_2017_02_rar_CCA_bray +theme_classic(base_size=15)+
  scale_color_viridis_d(name="Days after\ninjection",option="magma")+
  scale_shape_discrete(name="Well type")+xlim(-1.5,2)

####################################################################
#all together
#need to change some aes issues on initial plots


fig3_a$layers[[1]]$aes_params$size<-2
fig3_b$layers[[1]]$aes_params$size<-2

sample_data(ps.EVO_all_02_rar)$EVOYear<-as.factor(sample_data(ps.EVO_all_02_rar)$EVOYear)

ordinate(ps.EVO_all_02_rar, method = "NMDS", distance = "bray", weighted=FALSE)->ord.EVO_all_02_rar_NMDS_bray

ord<-ord.EVO_all_02_rar_NMDS_bray
choices<-c(1,2)
scaling<-1
pt.size<-3
angle<-20
len<-0.5
arrow.col<-"black"
unit<-"cm"
df_ord <- vegan::scores(ord, display = "sites", choices = choices, 
                        scaling = scaling)
df_ord <- as.data.frame(df_ord)
axis.labels <- ggordiplots::ord_labels(ord)[choices]
#there's some code here depending on "groups", which may be specified in envfit
#a formula may also be specified

#changing to the names used for plotting
colnames(df_ord) <- c("x", "y")
axis.labels <- ggordiplots::ord_labels(ord)[choices]

xlab <- "NMDS1"
ylab <- "NMDS2"

merge(as.data.frame(sample_data(ps.EVO_all_02_rar)), df_ord, by=0)->df_ord_2
df_ord_2->df_ord_both_combined

envfit(ordinate(subset_samples(ps.EVO_all_02_rar, DaysAfter!=372),
                method='NMDS', distance='bray'), 
       env=sample_data(subset_samples(ps.EVO_all_02_rar, DaysAfter!=372))[,c(13,15:18), drop=FALSE], na.rm=TRUE,
       permutations=999)->envfit.EVO_all_02_rar_NMDS_bray

df_arrows <- as.data.frame(vegan::scores(envfit.EVO_all_02_rar_NMDS_bray, "vectors"))
mult <- scale_arrow(df_arrows, df_ord_2[, c("x", "y")])
df_arrows <- mult * df_arrows
df_arrows$var <- rownames(df_arrows)
df_arrows$p.val <- envfit.EVO_all_02_NMDS_bray$vectors$pvals
colnames(df_arrows) <- c("x", "y", "var", "p.val")
df_arrows <- df_arrows[df_arrows$p.val <= 0.05, ]
df_arrows->df_arrows_all

fig3_c<-ggplot(subset(df_ord_both_combined, EVOYear=='2009'), aes(x=x,y=y))+
  geom_point(aes(shape=well_type,color=DaysAfterfactor),size=2)+
  scale_color_viridis_d(option='viridis')+theme_classic(base_size=15)+
  scale_shape_manual(values = c(1,16))+xlab('NMDS1')+ylab('NMDS2')+
  geom_text(data = envfit_coords_2, aes(x = NMDS1, y = NMDS2), colour = "black", 
            label = row.names(envfit_coords_2))+
  geom_segment(data = df_arrows_all, aes(x = 0, 
                                         xend = x, y = 0, yend = y), 
               arrow = arrow(angle = 20, 
                             length = unit(0.5, 'cm')), color = 'black')+
  ggnewscale::new_scale('shape')+ ggnewscale::new_scale_color()+
  geom_point(data=subset(df_ord_2, EVOYear=='2017'),aes(shape=well_type,color=DaysAfterfactor),size=2,stroke=1.05)+
  scale_color_viridis_d(option='magma')+
  scale_shape_manual(values = c(2,17))+
  ggtitle('c')+
  theme(text = element_text( size = 12), 
        plot.title = element_text(hjust = 1),
        legend.position = 'none')

get_legend(fig3_c+theme(legend.position='right'))->fig3legend

fig3_b+xlim(-2,2)->fig3_b

Fig_3<-gridExtra::grid.arrange(
  fig3_a +
    theme_classic(base_size=15)+
    scale_color_viridis_d(name="Days after\ninjection",option="viridis")+
    ggtitle('a')+
    theme(text = element_text( size = 12), 
          plot.title = element_text(hjust = 1),
          legend.position = 'none')+
    scale_shape_manual(name="Well type",values = c(1,16), labels=c('2009 Control','2009 Monitoring')),
  
  fig3_b +
    theme_classic(base_size=15)+
    scale_color_viridis_d(name="Days after\ninjection",option="plasma")+
    ggtitle('b')+
    theme(text = element_text( size = 12), 
          plot.title = element_text(hjust = 1),
          legend.position = 'none')+
    scale_shape_manual(values = c(2,17),name="Well type", labels=c('2017 Control','2017 Monitoring')),
  
  fig3_c,
  
  
  nrow=2,ncol=2#,
  #layout_matrix=rbind(c(1,1,2,2), c(3,3,NA,NA))
)

export::graph2ppt(file='Figure_3',upscale=T)


ggsave(plot=Fig_3,
       filename='Figure_3.tiff',
       device='tiff',
       width=190,
       height=190,
       units='mm')

ggsave(plot=get_legend(fig3_a +
                         theme_classic(base_size=15)+
                         scale_color_viridis_d(name="Days after\ninjection, 2009",option="viridis")+
                         ggtitle('a')+
                         theme(text = element_text( size = 12), 
                               plot.title = element_text(hjust = 1),
                               legend.position = 'bottom')+
                         scale_shape_manual(name="Well type",values = c(1,16), labels=c('2009 Control','2009 Monitoring'))),
       filename = 'fig3legend_a.tiff',
       dpi=1000,
       width=240,
       height=60,
       units='mm')
ggsave(plot=get_legend(fig3_b +
                         theme_classic(base_size=15)+
                         scale_color_viridis_d(name="Days after\ninjection, 2017",option="plasma")+
                         ggtitle('b')+
                         theme(text = element_text( size = 12), 
                               plot.title = element_text(hjust = 1),
                               legend.position = 'bottom')+
                         scale_shape_manual(values = c(2,17),name="Well type", labels=c('2017 Control','2017 Monitoring'))),
       filename = 'fig3legend_a.tiff',
       dpi=1000,
       width=240,
       height=60,
       units='mm')


#statistical analysis of beta diversity:
#ANOSIM
EVO_sampdata_09$anosim_group<-rep(NA,nrow(EVO_sampdata_09))
EVO_sampdata_17$anosim_group<-rep(NA,nrow(EVO_sampdata_17))
EVO_sampdata_09[which(EVO_sampdata_09$well_type=='Control'|EVO_sampdata_09$DaysAfter==-28),'anosim_group']<-'ctrl'
EVO_sampdata_09[which(EVO_sampdata_09$well_type!='Control'& EVO_sampdata_09$DaysAfter!=-28),'anosim_group']<-'mon'
EVO_sampdata_17[which(EVO_sampdata_17$well_type=='Control'|EVO_sampdata_17$DaysAfter==-6),'anosim_group']<-'ctrl'
EVO_sampdata_17[which(EVO_sampdata_17$well_type!='Control'& EVO_sampdata_17$DaysAfter!=-6),'anosim_group']<-'mon'


bray.EVO_2017_02<-distance(ps.EVO_2017_02_rar,method='bray')
bray.EVO_2009<-distance(ps.EVO_2009_rar,method='bray')
summary(anosim(x=bray.EVO_2009, grouping = EVO_sampdata_09$anosim_group))
summary(anosim(x=bray.EVO_2017_02, grouping = EVO_sampdata_17$anosim_group))

#PERMANOVA
adonis2(bray.EVO_2009 ~ group,
        data=data.frame(sample_data(ps.EVO_2009_rar)))

adonis2(bray.EVO_2017_02 ~ group,
        data=data.frame(sample_data(ps.EVO_2017_02_rar)))

#pairwise adonis 
#(not enough samples to call anything significant, but enough to give comparisons between days)
distance(ps.EVO_all_02_rar, method="bray")->bray_all
length(levels(as.factor(as.data.frame(sample_data(ps.EVO_all_02_rar))$group)))

new_group_names<- c("092con"="2009 control","092mon-28"="2009 -28", "092mon4"="2009 4","092mon17"="2009 17",
                    "092mon31"="2009 31", "092mon80"="2009 80",  "092mon140"="2009 140", "092mon269"="2009 269",
                    "172con"="2017 control","172mon-6"="2017 -6","172mon1"="2017 1", "172mon8"="2017 8",
                    "172mon15"="2017 15","172mon22"="2017 22","172mon50"="2017 50", "172mon78"="2017 78",
                    "172mon106"="2017 106", "172mon134"="2017 134","172mon372"="2017 372")

new_group_names_2<-c("092con"="2009\ncontrol","092mon-28"="-28", "092mon4"="4","092mon17"="17",
                     "092mon31"="31", "092mon80"="80",  "092mon140"="140", "092mon269"="269",
                     "172con"="2017\ncontrol","172mon-6"="-6","172mon1"="1", "172mon8"="8",
                     "172mon15"="15","172mon22"="22","172mon50"="50", "172mon78"="78",
                     "172mon106"="106", "172mon134"="134","172mon372"="372")


EVO_all_pairwiseadonisresults_complexmodel<-pairwiseAdonis::pairwise.adonis2(bray_all~group/EVOYear,
                                                                             data=data.frame(sample_data(ps.EVO_all_02_rar)),
                                                                             strata = data.frame(sample_data(ps.EVO_all_02_rar))$EVOYear)
betadivheatmapdf_complexmodel<-data.frame(matrix(ncol=19, nrow = 19))
colnames(betadivheatmapdf_complexmodel)<-levels(as.factor(as.data.frame(sample_data(ps.EVO_all_02_rar))$group))
rownames(betadivheatmapdf_complexmodel)<-levels(as.factor(as.data.frame(sample_data(ps.EVO_all_02_rar))$group))

for (i in 2:length(EVO_all_pairwiseadonisresults_complexmodel)){
  EVO_all_pairwiseadonisresults_complexmodel[[i]][1,3]->r2
  sub("_vs_.*",'',names(EVO_all_pairwiseadonisresults_complexmodel)[i] )->x
  sub(".*_vs_",'',names(EVO_all_pairwiseadonisresults_complexmodel)[i] )->y
  r2->betadivheatmapdf_complexmodel[x,y]
  r2->betadivheatmapdf_complexmodel[y,x]
}



levels(as.factor(sample_data(ps.EVO_all_02_rar)$group))->grouporder
grouporder[c(1,2,7,4,6,8,3,5,9,10, 11,19,14,15,17,18,12,13,16)]->grouporder
betadivheatmapdf_complexmodel[grouporder,grouporder]->betadivheatmapdf_reorder_complexmodel
pheatmap::pheatmap(betadivheatmapdf_reorder_complexmodel,
                   cluster_cols = F, cluster_rows=F ,na_col="white", border_color=NA)

dist.3col(as.dist(betadivheatmapdf_reorder_complexmodel))->betadivheatmapdf_plotting_complexmodel

betadivheatmapdf_plotting_complexmodel$name1<-factor(betadivheatmapdf_plotting_complexmodel$name1, levels=grouporder)
betadivheatmapdf_plotting_complexmodel$name2<-factor(betadivheatmapdf_plotting_complexmodel$name2, levels=grouporder)

#plotting
fig3_d<-ggplot(data=betadivheatmapdf_plotting_complexmodel, aes(name1, name2, fill=dis, size=dis,color=dis))+
  geom_point(stroke=0)+theme_classic(base_size = 15)+
  scale_color_viridis_c(limits=c(0.1, 0.7),
                        breaks=seq(0.1, 0.7, by=0.1), 
                        name=expression("PERMANOVA R"^"2"),
                        option='G')+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0))+
  scale_size_continuous(range=c(0.1,5),limits=c(0.1, 0.7), breaks=seq(0.1, 0.7, by=0.1), name=expression("PERMANOVA R"^"2"))+
  scale_color_continuous(limits=c(0.1, 0.7), breaks=seq(0.1, 0.7, by=0.1), name=expression("PERMANOVA R"^"2"))+
  scale_fill_continuous(limits=c(0.1, 0.7), breaks=seq(0.1, 0.7, by=0.1), name=expression("PERMANOVA R"^"2"))+
  guides(color= guide_legend(), size=guide_legend(), fill=guide_legend())+ 
  scale_x_discrete(labels=new_group_names_2)+ scale_y_discrete(labels=new_group_names_2)


#full figure
Fig_3<-gridExtra::grid.arrange(
  fig3_a +
    theme_classic(base_size=15)+
    scale_color_viridis_d(name="Days after\ninjection",option="viridis")+
    ggtitle('a')+
    theme(text = element_text( size = 12), 
          plot.title = element_text(hjust = 1),
          legend.position = 'none')+
    scale_shape_manual(name="Well type",values = c(1,16), labels=c('2009 Control','2009 Monitoring')),
  
  fig3_b +
    theme_classic(base_size=15)+
    scale_color_viridis_d(name="Days after\ninjection",option="plasma")+
    ggtitle('b')+
    theme(text = element_text( size = 12), 
          plot.title = element_text(hjust = 1),
          legend.position = 'none')+
    scale_shape_manual(values = c(2,17),name="Well type", labels=c('2017 Control','2017 Monitoring')),
  
  fig3_c,
  fig3_d+
    ggtitle('d')+
    theme(text = element_text( size = 12),
               legend.position='bottom', 
               plot.title = element_text(hjust = 1)),
  
  nrow=2,ncol=2#,
  #layout_matrix=rbind(c(1,1,2,2), c(3,3,NA,NA))
)

export::graph2ppt(file='Figure_3',upscale=T)


ggsave(plot=Fig_3,
       filename='Figure_3_withd.tiff',
       device='tiff',
       width=190,
       height=190,
       units='mm')

################################################################################
####################Figure S3: community average 16S rrn copy numbers
################################################################################
library(stringr)
#for the following code we do need to transform these to relative abundance
transform_sample_counts(ps.EVO_all_02_rar, function(x) x/sum(x))->ps.EVO_all_02_rar_RA
subset_samples(ps.EVO_all_02_rar_RA, EVOYear==2009)->ps.EVO_2009_rar_RA
subset_samples(ps.EVO_all_02_rar_RA, EVOYear==2017)->ps.EVO_2017_02_rar_RA
ghostbust(ps.EVO_2009_rar_RA)->ps.EVO_2009_rar_RA
ghostbust(ps.EVO_2017_02_rar_RA)->ps.EVO_2017_02_rar_RA


#first, import the rrn database from https://rrndb.umms.med.umich.edu/
rrndb <- read.delim("~/R_files/EVO_09vs17_inputs_usearch/for_data_sharing/rrnDB-5.8.tsv")
rrndb_reduced<-rrndb[,c(1:10,12),drop=F]
#and also new_phyla_names, which contains the name changes from the past few years
new_phyla_names <- read.delim("~/R_files/EVO_09vs17_inputs_usearch/for_data_sharing/new_phyla_names.txt")

#get the tax tables
data.frame(tax_table(ps.EVO_2009_rar))->EVO09tax
EVO09tax[,-c(8:13)]->EVO09tax

data.frame(tax_table(ps.EVO_2017_02_rar))->EVO17tax
EVO17tax[,-c(8:13)]->EVO17tax

EVO09tax$rrndb_copyno<-c(rep(NA, nrow(EVO09tax)))
EVO17tax$rrndb_copyno<-c(rep(NA, nrow(EVO17tax)))

#generate phylogenetic distance objects
cophenetic.phylo(phy_tree(ps.EVO_2009_rar))->phydist.EVO09_all
cophenetic.phylo(phy_tree(ps.EVO_2017_02_rar))->phydist.EVO17_all

#run the code to find 16S copy number per taxon
for (i in 1:nrow(EVO09tax)){
  
  #get the tax table entry for the ASV
  EVO09tax[i,,drop=FALSE]->temptax
  
  
  #this is code that checks if the ASV is only classified at the kingdom level
  #or if it's completely unclassified
  #then finds the nearest neighbor that is classified at at least the phylum level
  #and substitutes it for the unclassified one
  
  #using the cophenetic.phylo command to create a distance matrix of each ASV
  #then select the ASV and its distance values to every other ASV
  phydist.EVO09_all[,rownames(EVO09tax)[i],drop=FALSE]->phydistcompare
  j<-1
  
  #if the phylum to species is unclassified, we change the temptax entry to the nearest neighbor
  #this keeps going until we get to the jth neighbor that has a phylum or more
  while(all(is.na(temptax[,2:7]))){
    rownames(phydistcompare)[which(rank(phydistcompare[,1])==j)]->neighbor
    EVO09tax[neighbor,,drop=FALSE]->temptax
    j<-j+1
  }
  
  #find the finest taxrank of the nearest neighbor classified at the phylum level
  
  if(any(is.na(temptax[,2:7]))){
    colnames(temptax)[min(which(is.na(temptax)))-1]->taxrank
    temptax[, taxrank]->lowesttax
  } else if(any(is.na(temptax[,2:7]))==FALSE){
    taxrank<-"Species"
    lowesttax<-temptax[,taxrank]
  }
  taxrank_num<-which(colnames(temptax)==taxrank)
  if(taxrank_num==2){
    newphyrow<-which(apply(X=new_phyla_names,MARGIN = 1, FUN = function(X) str_which(X, lowesttax))==1)
    new_phyla_names[newphyrow,1:2]->lowesttax
  }
  #get the character vector of the lowest tax lvl 
  #this also makes the lowesttax into a multiple character
  as.character(lowesttax)->lowesttax
  
  #make lowesttax with multiple entries compatible with stringr via the following
  if(length(lowesttax)>1){
    paste(lowesttax, collapse = "|")->lowesttax
  }
  
  #check for matches in each of the areas it could be in the database
  stringr::str_which(rrndb_reduced$RDP.taxonomic.lineage, lowesttax)->RDP_taxno
  stringr::str_which(rrndb_reduced$RDP.taxa, lowesttax)->RDP_taxno_2
  stringr::str_which(rrndb_reduced$NCBI.scientific.name, lowesttax)->ncbi_taxno
  stringr::str_which(rrndb_reduced$Data.source.organism.name, lowesttax)->source_taxno
  
  #make a list of the vectors
  list(RDP_taxno, RDP_taxno_2, ncbi_taxno,source_taxno)->rrndb_taxnos_list
  
  #initialize the while loop by having n set to 2 so we can go up one taxonomic level
  taxrank_num<-which(colnames(temptax)==taxrank)
  
  #this loop checks to see if we couldn't find any matches to that taxonomic rank
  #then it goes up until it finds at least one match (ie the while loop closes)
  while(any(lengths(rrndb_taxnos_list)>0)==FALSE &&
        taxrank_num>=2){
    temptax[,taxrank_num]->lowesttax
    as.character(lowesttax)->lowesttax
    stringr::str_which(rrndb_reduced$RDP.taxonomic.lineage, lowesttax)->RDP_taxno
    stringr::str_which(rrndb_reduced$RDP.taxa, lowesttax)->RDP_taxno_2
    stringr::str_which(rrndb_reduced$NCBI.scientific.name, lowesttax)->ncbi_taxno
    stringr::str_which(rrndb_reduced$Data.source.organism.name, lowesttax)->source_taxno
    list(RDP_taxno, RDP_taxno_2, ncbi_taxno,source_taxno)->rrndb_taxnos_list
    
    taxrank_num<-taxrank_num-1
  }
  
  #now that we should have at least one match, get the row numbers in the rrndb
  unique(unlist(rrndb_taxnos_list[which(lengths(rrndb_taxnos_list)>0)]))->rrndb_rows
  
  #if that's still nothing, rrndb_rows should be NULL, so set up an if statement
  #and iterate if it's null
  if(is.null(rrndb_rows)==FALSE){
    mean(rrndb_reduced[rrndb_rows,"X16S.gene.count"])->mean_16sgeneabund
    mean_16sgeneabund->EVO09tax[i,"rrndb_copyno"]
  } else {
    next
  }
}


#adding the rrn to sample data
EVO09sampdata<-data.frame(sample_data(ps.EVO_2009_rar))
EVO09sampdata$rrncopyno<-rep(NA, nrow(EVO09sampdata))
which(rownames(EVO09tax) %in%
        rownames(otu_table(ps.EVO_2009_rar_RA)))->otureduce

#then start a for loop for every sample
for(i in 1:nsamples(ps.EVO_2009_rar_RA)){
  
  otu_table(ps.EVO_2009_rar_RA)[otureduce,i,drop=FALSE]->sampotus
  EVO09tax[-which(is.na(EVO09tax$rrndb_copyno)), ]->tax_reduced
  #only using taxa with values generated above
  tophalf<-sum(sampotus[which(rownames(sampotus)%in%rownames(tax_reduced))])
  bottomhalf<-c(0)
  for (j in 1:nrow(sampotus)){
    relabund<-sampotus[j,]
    rownames(sampotus)[j]->otuname
    copyno<-EVO09tax[otuname,"rrndb_copyno"]
    if(is.na(copyno)==TRUE){next}
    as.numeric(relabund)/as.numeric(copyno)->sovern
    bottomhalf<-bottomhalf+sovern
  }
  as.numeric(tophalf)/as.numeric(bottomhalf)->EVO09sampdata[colnames(sampotus),"rrncopyno"]
}

#plotting
EVO09sampdata[,c(1:5,19,22)]->EVO09rrncopydata
p.EVO09_16Scopynumber<-ggplot(data=EVO09rrncopydata, aes(x=DaysAfterfactor, y=rrncopyno))+ geom_boxplot()+geom_point(aes(color=well_type))+
  theme_classic(base_size=15)+xlab("Days after injection")+ylab("rrn")+labs(color="Well type")

p.EVO09_16Scopynumber
##############################
#now 2017
for (i in 1:nrow(EVO17tax)){
  
  #get the tax table entry for the ASV
  EVO17tax[i,,drop=FALSE]->temptax
  
  
  #this is code that checks if the ASV is only classified at the kingdom level
  #or if it's completely unclassified
  #then finds the nearest neighbor that is classified at at least the phylum level
  #and substitutes it for the unclassified one
  
  #using the cophenetic.phylo command to create a distance matrix of each ASV
  #then select the ASV and its distance values to every other ASV
  phydist.EVO17_all[,rownames(EVO17tax)[i],drop=FALSE]->phydistcompare
  j<-1
  #if the phylum to species is unclassified, we change the temptax entry to the nearest neighbor
  #this keeps going until we get to the jth neighbor that has a phylum or more
  while(all(is.na(temptax[,2:7]))){
    rownames(phydistcompare)[which(rank(phydistcompare[,1])==j)]->neighbor
    EVO17tax[neighbor,,drop=FALSE]->temptax
    j<-j+1
  }
  
  #find the finest taxrank of the nearest neighbor classified at the phylum level
  
  if(any(is.na(temptax[,2:7]))){
    colnames(temptax)[min(which(is.na(temptax)))-1]->taxrank
    temptax[, taxrank]->lowesttax
  } else if(any(is.na(temptax[,2:7]))==FALSE){
    taxrank<-"Species"
    lowesttax<-temptax[,taxrank]
  }
  taxrank_num<-which(colnames(temptax)==taxrank)
  if(taxrank_num==2){
    newphyrow<-which(apply(X=new_phyla_names,MARGIN = 1, FUN = function(X) str_which(X, lowesttax))==1)
    new_phyla_names[newphyrow,1:2]->lowesttax
  }
  #get the character vector of the lowest tax lvl 
  #this also makes the lowesttax into a multiple character
  as.character(lowesttax)->lowesttax
  
  #make lowesttax with multiple  compatible with stringr via the following
  if(length(lowesttax)>1){
    paste(lowesttax, collapse = "|")->lowesttax
  }
  
  #check for matches in each of the areas it could be in the database
  stringr::str_which(rrndb_reduced$RDP.taxonomic.lineage, lowesttax)->RDP_taxno
  stringr::str_which(rrndb_reduced$RDP.taxa, lowesttax)->RDP_taxno_2
  stringr::str_which(rrndb_reduced$NCBI.scientific.name, lowesttax)->ncbi_taxno
  stringr::str_which(rrndb_reduced$Data.source.organism.name, lowesttax)->source_taxno
  
  #make a list of the vectors
  list(RDP_taxno, RDP_taxno_2, ncbi_taxno,source_taxno)->rrndb_taxnos_list
  
  #initialize the while loop by having n set to 2 so we can go up one taxonomic level
  taxrank_num<-which(colnames(temptax)==taxrank)
  #this loop checks to see if we couldn't find any matches to that taxonomic rank
  #then it goes up until it finds at least one match (ie the while loop closes)
  while(any(lengths(rrndb_taxnos_list)>0)==FALSE &&
        taxrank_num>=2){
    temptax[,taxrank_num]->lowesttax
    as.character(lowesttax)->lowesttax
    stringr::str_which(rrndb_reduced$RDP.taxonomic.lineage, lowesttax)->RDP_taxno
    stringr::str_which(rrndb_reduced$RDP.taxa, lowesttax)->RDP_taxno_2
    stringr::str_which(rrndb_reduced$NCBI.scientific.name, lowesttax)->ncbi_taxno
    stringr::str_which(rrndb_reduced$Data.source.organism.name, lowesttax)->source_taxno
    list(RDP_taxno, RDP_taxno_2, ncbi_taxno,source_taxno)->rrndb_taxnos_list
    
    taxrank_num<-taxrank_num-1
  }
  
  #now that we should have at least one match, get the row numbers in the rrndb
  unique(unlist(rrndb_taxnos_list[which(lengths(rrndb_taxnos_list)>0)]))->rrndb_rows
  
  #if that's still nothing, rrndb_rows should be NULL, so set up an if statement
  #and iterate if it's null
  if(is.null(rrndb_rows)==FALSE){
    mean(rrndb_reduced[rrndb_rows,"X16S.gene.count"])->mean_16sgeneabund
    mean_16sgeneabund->EVO17tax[i,"rrndb_copyno"]
  } else {
    next
  }
}


#make sure taxa_are_rows==TRUE
taxa_are_rows(ps.EVO_2017_02_rar_RA)

EVO17sampdata<-as.data.frame(sample_data(ps.EVO_2017_02_rar))
EVO17sampdata$rrncopyno<-rep(NA, nrow(EVO17sampdata))
which(rownames(EVO17tax) %in%
        rownames(otu_table(ps.EVO_2017_02_rar_RA)))->otureduce

#then start a for loop for every sample
as.data.frame(otu_table(ps.EVO_2017_02_rar_RA)[otureduce,,drop=FALSE])->sampotus_forchecking

for(i in 1:nsamples(ps.EVO_2017_02_rar_RA)){
  otu_table(ps.EVO_2017_02_rar_RA)[otureduce,i,drop=FALSE]->sampotus
  EVO17tax[-which(is.na(EVO17tax$rrndb_copyno)), ]->tax_reduced
  tophalf<-sum(sampotus[which(rownames(sampotus)%in%rownames(tax_reduced))])
  bottomhalf<-c(0)
  for (j in 1:nrow(sampotus)){
    relabund<-sampotus[j,]
    rownames(sampotus)[j]->otuname
    copyno<-EVO17tax[otuname,"rrndb_copyno"]
    if(is.na(copyno)==TRUE){next}
    as.numeric(relabund)/as.numeric(copyno)->sovern
    bottomhalf<-bottomhalf+sovern
  }
  as.numeric(tophalf)/as.numeric(bottomhalf)->EVO17sampdata[colnames(sampotus),"rrncopyno"]
}

#plotting
EVO17sampdata[,c(1:5,19,22)]->EVO17rrncopydata
p.EVO17_16Scopynumber<-ggplot(data=EVO17rrncopydata, aes(x=DaysAfterfactor, y=rrncopyno))+ geom_boxplot()+geom_point(aes(color=well_type))+
  theme_classic(base_size=15)+xlab("Days after injection")+ylab("rrn")+labs(color="Well type")
p.EVO17_16Scopynumber


rbind(EVO09rrncopydata,EVO17rrncopydata)->EVOall_rrn

annotation_text<-data.frame(DaysAfterfactor=factor('372',levels=levels(EVOall_rrn$DaysAfterfactor)), rrncopyno=7, lab='a', EVOYear=factor('2009', levels=c('2009','2017')))

Fig_s3<-ggplot(data=EVOall_rrn, aes(x=DaysAfterfactor, y=rrncopyno))+ geom_boxplot()+geom_point(aes(color=well_type))+
  theme_classic(base_size=15)+xlab("Days after injection")+ylab("rrn")+labs(color="Well type")+
  facet_wrap(~EVOYear, scales = 'free_x')+theme(strip.background = element_blank(),
                                                strip.text=element_blank(),
                                                legend.position='bottom')
Fig_s3

export::graph2ppt(upscale=T)
ggsave(filename='Figure_S3.tiff',
       plot=)


#################################################################################
#################################################################################
####################Figure 4 & Table S2 (Differential abundance analysis)
#################################################################################
#################################################################################
library(Maaslin2)
library(treemapify)


#running MaAslin2 with mostly defaults
mal.EVO_2009<-Maaslin2(input_data = data.frame(otu_table(ps.EVO_2009_rar)),
                       input_metadata = data.frame(sample_data(ps.EVO_2009_rar)),
                       output = "~/R_files/EVO_09vs17_outputs_usearch", #change for your own data
                       min_abundance = 0.01,  #relative abundance, 1% cutoff
                       min_prevalence = 0.125,    #present in 1 timepoint
                       normalization = "TSS", #relative abundance  
                       transform = "LOG",
                       analysis_method = "LM",
                       max_significance = 0.25,
                       fixed_effects = "DaysAfterfactor",
                       random_effects = "WellID",
                       correction = "BH",
                       standardize = FALSE,
                       plot_scatter = FALSE,
                       plot_heatmap = FALSE,
                       cores = 8,
                       reference = "DaysAfterfactor,-28")


mal.EVO_2017<-Maaslin2(input_data = data.frame(otu_table(ps.EVO_2017_02_rar)),
                       input_metadata = data.frame(sample_data(ps.EVO_2017_02_rar)),
                       output = "~/R_files/EVO_09vs17_outputs_usearch", 
                       min_abundance = 0.01,  #relative abundance, 1% 
                       min_prevalence = 0.1,    #pretty strict, present in 2 timepoints
                       normalization = "TSS", #default 
                       transform = "LOG", #default when used with LM
                       analysis_method = "LM",  #default
                       max_significance = 0.25, #default
                       fixed_effects = "DaysAfterfactor",
                       random_effects = "WellID", #using wellID as a random effect
                       correction = "BH",
                       standardize = FALSE,
                       plot_scatter = FALSE,
                       plot_heatmap = FALSE,
                       cores = 8,
                       reference = "DaysAfterfactor,-6")

#identify significantly DA features at 0.05 level
DAstimfeat09<-unique(mal.EVO_2009$results$feature[which(mal.EVO_2009$results$qval<0.05)])
DAstimfeat17<-unique(mal.EVO_2017$results$feature[which(mal.EVO_2017$results$qval<0.05)])

#maaslin2 adds an x to feature names that start with a number, so we have to remove them
#making a function for convenience
removeX<-function(asvnames){
  for (i in 1:length(asvnames)){
    if(nchar(asvnames)[i]==33){
      substring(asvnames[i],2)->asvnames[i]
    }
    
  }
  return(asvnames)
}

removeX(DAstimfeat09)->DAstimfeat09
removeX(DAstimfeat17)->DAstimfeat17

#make the results wide then correct the feature names
mal.EVO_2009$results->mal.EVO_2009_res
mal.EVO_2009_res[,c("feature",'value','coef','stderr','qval')]->mal.EVO_2009_res_reduced
DAstimres09<-as.data.frame(pivot_wider(mal.EVO_2009_res_reduced, names_from = 2, values_from =3:5))
removeX(DAstimres09$feature)->DAstimres09$feature

mal.EVO_2017$results->mal.EVO_2017_res
mal.EVO_2017_res[,c("feature",'value','coef','stderr','qval')]->mal.EVO_2017_res_reduced
DAstimres17<-as.data.frame(pivot_wider(mal.EVO_2017_res_reduced, names_from = 2, values_from =3:5))
removeX(DAstimres17$feature)->DAstimres17$feature

#see how many ASVs overlap in Maaslin2 results
length(intersect(DAstimres09$feature,DAstimres17$feature))

#reduce to those with qval<0.05 at one time point or another
rownames(DAstimres09)<-DAstimres09$feature
rownames(DAstimres17)<-DAstimres17$feature

DAstimres09[DAstimfeat09,]->DAstimres09_sig
DAstimres17[DAstimfeat17,,drop=F]->DAstimres17_sig
length(intersect(DAstimres09_sig$feature,DAstimres17_sig$feature))

#assess those in phase 2, the reductive phase
#using a coeff cutoff of 1, which is an e-fold change (2.7x)
#I've also tried a 0 cutoff, which is a 1-fold change or more (ie not decreased from preinject)
DAstimres09_sig_4173180<-DAstimres09_sig[which(DAstimres09_sig[,'coef_4']>1|
                                                 DAstimres09_sig[,'coef_17']>1|
                                                 DAstimres09_sig[,'coef_31']>1|
                                                 DAstimres09_sig[,'coef_80']>1),,drop=F]
#could also omit day 80
DAstimres09_sig_41731<-DAstimres09_sig[which(DAstimres09_sig[,'coef_4']>1|
                                               DAstimres09_sig[,'coef_17']>1|
                                               DAstimres09_sig[,'coef_31']>1),,drop=F]

DAstimres17_sig_8152250<-DAstimres17_sig[which(DAstimres17_sig[,'coef_8']>1|
                                                 DAstimres17_sig[,'coef_15']>1|
                                                 DAstimres17_sig[,'coef_22']>1|
                                                 DAstimres17_sig[,'coef_50']>1),,drop=F]

nrow(DAstimres09_sig_4173180)
nrow(DAstimres17_sig_8152250)
length(intersect(rownames(DAstimres17_sig_8152250), rownames(DAstimres09_sig_4173180)))

merge(as.data.frame(tax_table(ps.EVO_2009_rar)[rownames(DAstimres09_sig_4173180),,drop=F]),
      DAstimres09_sig_4173180, by=0)->DAstimactive09
merge(as.data.frame(tax_table(ps.EVO_2017_02_rar)[rownames(DAstimres17_sig_8152250),,drop=F]),
      DAstimres17_sig_8152250, by=0)->DAstimactive17
rownames(DAstimactive09)<-DAstimactive09$Row.names
rownames(DAstimactive17)<-DAstimactive17$Row.names
DAstimactive09[,-1]->DAstimactive09
DAstimactive17[,-1]->DAstimactive17

#subset phyloseq objects for further analysis
ps.EVO_2009_rar_RA_stimactive<-subset_taxa(ps.EVO_2009_rar_RA, 
                                           taxa_names(ps.EVO_2009_rar_RA) %in% rownames(DAstimactive09))

ps.EVO_2017_02_rar_RA_stimactive<-subset_taxa(ps.EVO_2017_02_rar_RA, 
                                              taxa_names(ps.EVO_2017_02_rar_RA) %in% rownames(DAstimactive17))

apply(as.data.frame(otu_table(ps.EVO_2009_rar_RA_stimactive)), 1, max)->maxRA09stimactive
merge(DAstimactive09, as.data.frame(maxRA09stimactive), by=0)->DAstimactive09_fortreemap
DAstimactive09_fortreemap$Row.names->rownames(DAstimactive09_fortreemap)
DAstimactive09_fortreemap[,-1]->DAstimactive09_fortreemap
DAstimactive09_fortreemap$maxcoef<-apply(DAstimactive09_fortreemap[,c(16:18,20)],1,max)
DAstimactive09_fortreemap[is.na(DAstimactive09_fortreemap)] <- "Unclassified"


################################################################################
#ITOL export
################################################################################
#these dataframes were initially used for treemap representations of the data
#instead, they will be exported to itol to construct figures using itol.toolkit

#prepare for export 
library(itol.toolkit) # main package
library(dplyr) # data manipulation
library(data.table) # file read
library(ape) # tree operation
library(stringr) # string operation
library(tidyr) # data manipulation

asvforitol<-unique(union(rownames(DAstimactive09_fortreemap),rownames(DAstimactive17_fortreemap)))

ps.EVO_all_02_rar_foritol<-subset_taxa(ps.EVO_all_02_rar, taxa_names(ps.EVO_all_02_rar) %in% asvforitol)
ape::write.tree(phy_tree(ps.EVO_all_02_rar_foritol),file='EVO_tree_foritol.nwk')

taxtabdf_foritol_EVO<-data.frame(tax_table(ps.EVO_all_02_rar_foritol))

#ggtree testing; having compatibility issues with it and/or tidytree, so I won't be using it
#leaving the code since some was used downstream

colnames(DAstimactive09_fortreemap)[34:35]<-c('maxcoef_09','Genus_overlap_09')
DAstimactive09_fortreemap[,c(1:7,33:35)]->EVO09treedat
EVO09treedat$ID<-rownames(EVO09treedat)
EVO09treedat<-EVO09treedat[,c(ncol(EVO09treedat),1:ncol(EVO09treedat)-1)]

colnames(DAstimactive17_fortreemap)[43:44]<-c('maxcoef_17','Genus_overlap_17')
DAstimactive17_fortreemap[,c(1:7,42:44)]->EVO17treedat
EVO17treedat$ID<-rownames(EVO17treedat)
EVO17treedat<-EVO17treedat[,c(ncol(EVO17treedat),1:ncol(EVO17treedat)-1)]

EVOalltreedat<-merge(EVO09treedat,EVO17treedat,all = T)
rownames(EVOalltreedat)<-EVOalltreedat$ID

taxtabEVOforitol<-data.frame(tax_table(ps.EVO_all_02_rar_foritol))
taxtabEVOforitol$ID<-rownames(taxtabEVOforitol)
taxtabEVOforitol[,c('ID','Class')]

grp_EVOtree<-list()
for(i in 1:length(unique(taxtabEVOforitol$Class))){
  grp_EVOtree[[length(grp_EVOtree)+1]]<-taxtabEVOforitol[which(taxtabEVOforitol$Class==unique(taxtabEVOforitol$Class)[i]),'ID']
  names(grp_EVOtree)[i]<-unique(taxtabEVOforitol$Class)[i]
}


EVOtree_DAstim<-ggtree(phy_tree(ps.EVO_all_02_rar_foritol), layout = 'circular', branch.length = 'none')

EVOtree_DAstim+groupOTU(EVOtree_DAstim,grp_EVOtree ,'Class')

EVOtree_DAstim%<+%
  EVOalltreedat +
  geom_tippoint(aes(x=x+1,size=maxRA09stimactive, color=maxcoef_09), shape=16, na.rm=T)+
  scale_colour_viridis_c(option='plasma')+
  ggnewscale::new_scale_colour()+
  geom_tippoint(aes(x=x+3,size=maxRA17stimactive, color=maxcoef_17), shape=17, na.rm=T)+
  scale_colour_viridis_c(option='viridis')


#itol code 

hub_EVODAstim<-create_hub(tree=phy_tree(ps.EVO_all_02_rar_foritol))

EVOalltreedat$logmaxRA09<-log(EVOalltreedat$maxRA09stimactive)
EVOalltreedat$logmaxRA17<-log(EVOalltreedat$maxRA17stimactive)
EVOalltreedat_noNA<-EVOalltreedat
EVOalltreedat_noNA[is.na(EVOalltreedat_noNA)]<-'X'

#cutoff: max RA>0.005

ps.EVO_all_02_rar_foritol_2<-subset_taxa(ps.EVO_all_02_rar_foritol, 
                                         taxa_names(ps.EVO_all_02_rar_foritol) %in% 
                                           rownames(EVOalltreedat)[which(EVOalltreedat$maxRA17stimactive>0.005 | EVOalltreedat$maxRA09stimactive>0.005)])
ape::write.tree(phy_tree(ps.EVO_all_02_rar_foritol_2),file='EVO_tree_foritol.nwk')

hub_EVODAstim<-create_hub(tree=phy_tree(ps.EVO_all_02_rar_foritol_2))


EVOalltreedat_noNA<-EVOalltreedat[which(EVOalltreedat$maxRA17stimactive>0.005 | EVOalltreedat$maxRA09stimactive>0.005),]
EVOalltreedat_noNA[is.na(EVOalltreedat_noNA)]<-'X'



unit_phylum_EVO<-create_unit(data = EVOalltreedat_noNA %>% select(ID, Phylum), 
                             key = "phylum_treecols", 
                             type = "TREE_COLORS", 
                             subtype = "range", 
                             tree=phy_tree(ps.EVO_all_02_rar_foritol_2))
unit_hmp_coefs <- create_unit(data = EVOalltreedat_noNA %>% 
                                select(ID, maxcoef_09,maxcoef_17), 
                              key = "coefs_hmp",
                              type = "DATASET_HEATMAP", 
                              tree=phy_tree(ps.EVO_all_02_rar_foritol_2))

unit_hmp_maxRA <- create_unit(data = EVOalltreedat_noNA %>% 
                                select(ID,logmaxRA09,logmaxRA17), 
                              key = "maxRA_hmp",
                              type = "DATASET_HEATMAP", 
                              tree=phy_tree(ps.EVO_all_02_rar_foritol_2))



hub_EVODAstim<-hub_EVODAstim+
  unit_phylum_EVO+
  unit_hmp_coefs+
  unit_hmp_maxRA

write_hub(hub_EVODAstim)

unit_Class_EVO<-create_unit(data = EVOalltreedat_noNA %>% select(ID, Class), 
                            key = "class_treecols", 
                            type = "TREE_COLORS", 
                            subtype = "range", 
                            tree=phy_tree(ps.EVO_all_02_rar_foritol_2))
write_unit(unit_Class_EVO)

unit_ordlab_EVO<-create_unit(data = EVOalltreedat_noNA %>% select(ID, Order), 
                             key = "order_labels", 
                             type = "LABELS", 
                             subtype = "range", 
                             tree=phy_tree(ps.EVO_all_02_rar_foritol_2))
write_unit(unit_ordlab_EVO)

unit_famlab_EVO<-create_unit(data = EVOalltreedat_noNA %>% select(ID, Family), 
                             key = "family_labels", 
                             type = "LABELS", 
                             subtype = "range", 
                             tree=phy_tree(ps.EVO_all_02_rar_foritol_2))
write_unit(unit_famlab_EVO)

#highlighting specific taxa from the paper with tree color brackets
#taxa: genus pelosinus, genus NK4A214_group, genus Desulforegula, genus Desulfovibrio, order Geobacterales, class Gammaproteobacteria
EVOalltreedat_noNA$Highlight<-rep(NA,nrow(EVOalltreedat_noNA))
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Genus=='Pelosinus'),'Highlight']<-'Pelosinus'
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Genus=='NK4A214_group'),'Highlight']<-'NK4A214_group'
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Genus=='Desulforegula'),'Highlight']<-'Desulforegula'
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Genus=='Desulfovibrio'),'Highlight']<-'Desulfovibrio'
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Order=='Geobacterales'),'Highlight']<-'Geobacterales'
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Class=='Gammaproteobacteria'),'Highlight']<-'Gammaproteobacteria'
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Class=='Gammaproteobacteria'),'Highlight']<-'Gammaproteobacteria'

#need to do all oscillospiraceae for NK4A214
EVOalltreedat_noNA[which(EVOalltreedat_noNA$Order=='Oscillospirales'),'Highlight']<-'Oscillospirales'

unit_taxofint_treecols<-create_unit(data=EVOalltreedat_noNA %>% select(ID,Highlight) %>%
                                      filter(is.na(Highlight)==F),
                                    key='highlight_treecols',
                                    type='TREE_COLORS',
                                    subtype='range',
                                    tree=phy_tree(ps.EVO_all_02_rar_foritol_2))
write_unit(unit_taxofint_treecols)

###############################################################################
################################################################################
################################################################################
################################################################################
#Fig 5 and Fig S4 community assembly mechanisms
################################################################################
################################################################################
library(NST)
library(iCAMP)
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
################################################################################
#icamp 09
comm.icamp.09<-t(data.frame(otu_table(ps.EVO_2009_rar)))
treat.icamp.09<-data.frame(sample_data(ps.EVO_2009_rar))
env.icamp.09<-data.frame(sample_data(ps.EVO_2009_rar))[, 7:18,drop=FALSE]
tree.icamp.09<-phy_tree(ps.EVO_2009_rar)

#save.wd.09<-"~/R_files/EVO_09vs17_inputs_usearch/for_data_sharing/iCAMP_outputs_09"
prefix.09<-"EVO09"
rand.time<-1000

setwd(save.wd.09)

pd.big.09=iCAMP::pdist.big(tree = tree.icamp.09, wd=save.wd.09, nworker = 8)
# output files:
# path.rda: a R object to list all the nodes and  edge lengthes from root to every tip. saved in R data format. an intermediate output when claculating phylogenetic distance matrix.
# pd.bin: BIN file (backingfile) generated by function big.matrix in R package bigmemory. This is the big matrix storing pairwise phylogenetic distance values. By using this bigmemory format file, we will not need memory but hard disk when calling big matrix for calculation.
# pd.desc: the DESC file (descriptorfile) to hold the backingfile (pd.bin) description.
# pd.taxon.name.csv: comma delimited csv file storing the IDs of tree tips (OTUs), serving as the row/column names of the big phylogenetic distance matrix.



ds = 0.2 
bin.size.limit = 24 


icres.09_bin24<- icamp.big(comm=comm.icamp.09, pd.desc = pd.big.09$pd.file, pd.spname=pd.big.09$tip.label,
                           pd.wd = pd.big.09$pd.wd, rand = 1000, tree=tree.icamp.09,
                           prefix = prefix.09, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                           phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                           phylo.metric = "bMPD", sig.index="Confidence", bin.size.limit = bin.size.limit, 
                           nworker = 8, rtree.save = FALSE, detail.save = TRUE, 
                           qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd.09, 
                           correct.special = TRUE, unit.sum = rowSums(comm.icamp.09), special.method = "depend",
                           ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
save(icres.09, file = "icres_09_bin24.rda")
bin.size.limit<-12
icres.09_bin12<- icamp.big(comm=comm.icamp.09, pd.desc = pd.big.09$pd.file, pd.spname=pd.big.09$tip.label,
                           pd.wd = pd.big.09$pd.wd, rand = 1000, tree=tree.icamp.09,
                           prefix = "EVO09bin12", ds = 0.2, pd.cut = NA, sp.check = TRUE,
                           phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                           phylo.metric = "bMPD", sig.index="Confidence", bin.size.limit = bin.size.limit, 
                           nworker = 8, rtree.save = FALSE, detail.save = TRUE, 
                           qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd.09, 
                           correct.special = TRUE, unit.sum = rowSums(comm.icamp.09), special.method = "depend",
                           ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
bin.size.limit<-48
icres.09_bin48<- icamp.big(comm=comm.icamp.09, pd.desc = pd.big.09$pd.file, pd.spname=pd.big.09$tip.label,
                           pd.wd = pd.big.09$pd.wd, rand = 1000, tree=tree.icamp.09,
                           prefix = "EVO09bin48", ds = 0.2, pd.cut = NA, sp.check = TRUE,
                           phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                           phylo.metric = "bMPD", sig.index="Confidence", bin.size.limit = bin.size.limit, 
                           nworker = 8, rtree.save = FALSE, detail.save = TRUE, 
                           qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd.09, 
                           correct.special = TRUE, unit.sum = rowSums(comm.icamp.09), special.method = "depend",
                           ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)

treat.icamp.09.boot<-treat.icamp.09[,"group", drop=FALSE]

#checking icamp.boot
icres.09.boot_bin24<-iCAMP::icamp.boot(icamp.result = icres.09_bin24$CbMPDiCBraya,
                                       treat=treat.icamp.09.boot)

icres.09.boot_bin12<-iCAMP::icamp.boot(icamp.result = icres.09_bin12$CbMPDiCBraya,
                                       treat=treat.icamp.09.boot)
icres.09.boot_bin48<-iCAMP::icamp.boot(icamp.result = icres.09_bin48$CbMPDiCBraya,
                                       treat=treat.icamp.09.boot)

icresboottest<-icres.09.boot_bin48
as.numeric(icresboottest$summary$Observed)->icresboottest$summary$Observed

icresboottest$summary[,1:3]->icresboottestplot
icresboottestplot$Process<-as.factor(icresboottestplot$Process)
icresboottestplot$Process<-recode_factor(icresboottestplot$Process,
                                         'Dispersal.Limitation'='DL',
                                         'Drift.and.Others'='DR',
                                         'Homogenizing.Dispersal'='HD',
                                         "Heterogeneous.Selection"="HeS",
                                         'Homogeneous.Selection'='HoS')


#plotting

plot.icres.09.boot.bin48<-ggplot2::ggplot(data=icresboottestplot,
                                          aes(fill=Process,
                                              x=factor(Group,levels = c("092con",
                                                                        "092mon-28",
                                                                        '092mon4',
                                                                        '092mon17',
                                                                        '092mon31',
                                                                        '092mon80',
                                                                        '092mon140',
                                                                        '092mon169')),
                                              y=Observed))+
  geom_bar(position = "stack", stat="identity") +xlab("Days after injection") +
  theme_minimal() +ylab("Observed relative importance")

p.icres.09<-plot.icres.09.boot.bin24+scale_fill_manual(values=cbPalette[2:6])+
  scale_x_discrete(labels=c('Control','-28','4','17','31','80','140','269'))+
  theme(text = element_text(size = 20))
p.icres.09+ ggtitle("(a)")+theme(plot.title = element_text(hjust = 1.2, vjust=-3))
export::graph2ppt(file='icamp09')

ggsave(filename = "iCAMP09.jpg",
       plot=p.icres.09,
       path="~/R_files/ggsave_destination")
########
#now icamp cate to see the differences between stimulated nd non-stimulated taxa

cate.09<-data.frame(matrix(NA, nrow=ntaxa(ps.EVO_2009_rar), ncol = 2))
cate.09[,1]<-taxa_names(ps.EVO_2009_rar)
cate.09[which(cate.09[,1] %in% taxa_names(ps.EVO_2009_rar_RA_stimactive)), 2]<-"stimulated"
cate.09[which(is.na(cate.09[,2])),2]<-"not_stimulated"
rownames(cate.09)<-cate.09[,1]
cate.09[,-1,drop=FALSE]->cate.09


clas.09<-data.frame(tax_table(ps.EVO_2009_rar))
icresbins.EVO_09<-icamp.bins(icres.09_bin24$detail, treat = NULL,
                             clas =clas.09 )

icrescate.EVO_09<-icamp.cate(icresbins.EVO_09,comm = icres.09_bin24$detail$comm,
                             cate = cate.09, treat = treat.icamp.09.boot)

icrescate.EVO_09$Ptx->iccate.09
iccate.09[,-c(1:2,9,15)]->iccate.09
iccate.09<-data.frame(pivot_longer(data=iccate.09, cols=2:11, names_to = c("Stimulation","Process"),
                                   names_sep = "[.]"))

plot.iccate.09<-ggplot2::ggplot(data=iccate.09, 
                                aes(fill=Process,
                                    x=factor(Group,levels = c("092con",
                                                              "092mon-28",
                                                              '092mon4',
                                                              '092mon17',
                                                              '092mon31',
                                                              '092mon80',
                                                              '092mon140',
                                                              '092mon169')), 
                                    y=value)) +
  geom_bar(position = "stack", stat="identity")+xlab("Days after injection")+
  facet_grid(.~Stimulation, labeller= as_labeller(c('not_stimulated'="Not stimulated",'stimulated'='Stimulated')))+
  theme_minimal() +ylab("Observed relative importance")

p.icrescate.09<-plot.iccate.09+scale_fill_manual(values=cbPalette[2:6])+
  scale_x_discrete(labels=c('Control','-28','4','17','31','80','140','269'))+
  theme(text = element_text(size = 20))
p.icrescate.09+ ggtitle("(c)")+theme(plot.title = element_text(hjust = 1.2, vjust=-3),
                                     axis.text.x = element_text(size=14,angle=45,vjust = 1.5,  hjust=1.3))
export::graph2ppt(file='icamp_category09')

ggsave(filename = "iCAMPcate09.jpg",
       plot=p.icrescate.09,
       path="~/R_files/ggsave_destination")



################################################################################
################################################################################
#same procedure for 2017
#icamp 17


bin.size.limit<-24
icres.17_bin24<- icamp.big(comm=comm.icamp.17, pd.desc = pd.big.17$pd.file, pd.spname=pd.big.17$tip.label,
                           pd.wd = pd.big.17$pd.wd, rand = 1000, tree=tree.icamp.17,
                           prefix = prefix.17, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                           phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                           phylo.metric = "bMPD", sig.index="Confidence", bin.size.limit = bin.size.limit, 
                           nworker = 8, rtree.save = FALSE, detail.save = TRUE, 
                           qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd.17, 
                           correct.special = TRUE, unit.sum = rowSums(comm.icamp.17), special.method = "depend",
                           ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)

bin.size.limit<-12
icres.17_bin12<- icamp.big(comm=comm.icamp.17, pd.desc = pd.big.17$pd.file, pd.spname=pd.big.17$tip.label,
                           pd.wd = pd.big.17$pd.wd, rand = 1000, tree=tree.icamp.17,
                           prefix = "EVO17bin12", ds = 0.2, pd.cut = NA, sp.check = TRUE,
                           phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                           phylo.metric = "bMPD", sig.index="Confidence", bin.size.limit = bin.size.limit, 
                           nworker = 8, rtree.save = FALSE, detail.save = TRUE, 
                           qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd.17, 
                           correct.special = TRUE, unit.sum = rowSums(comm.icamp.17), special.method = "depend",
                           ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
bin.size.limit<-48
icres.17_bin48<- icamp.big(comm=comm.icamp.17, pd.desc = pd.big.17$pd.file, pd.spname=pd.big.17$tip.label,
                           pd.wd = pd.big.17$pd.wd, rand = 1000, tree=tree.icamp.17,
                           prefix = "EVO17bin48", ds = 0.2, pd.cut = NA, sp.check = TRUE,
                           phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                           phylo.metric = "bMPD", sig.index="Confidence", bin.size.limit = bin.size.limit, 
                           nworker = 8, rtree.save = FALSE, detail.save = TRUE, 
                           qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd.17, 
                           correct.special = TRUE, unit.sum = rowSums(comm.icamp.17), special.method = "depend",
                           ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)

treat.icamp.17.boot<-treat.icamp.17[,"group", drop=FALSE]

icres.17.boot_bin24<-iCAMP::icamp.boot(icamp.result = icres.17_bin24$CbMPDiCBraya,
                                       treat=treat.icamp.17.boot)

icres.17.boot_bin12<-iCAMP::icamp.boot(icamp.result = icres.17_bin12$CbMPDiCBraya,
                                       treat=treat.icamp.17.boot)
icres.17.boot_bin48<-iCAMP::icamp.boot(icamp.result = icres.17_bin48$CbMPDiCBraya,
                                       treat=treat.icamp.17.boot)




icresboottest<-icres.17.boot_bin24
as.numeric(icresboottest$summary$Observed)->icresboottest$summary$Observed

icresboottest$summary[,1:3]->icresboottestplot
icresboottestplot$Process<-as.factor(icresboottestplot$Process)
icresboottestplot$Process<-recode_factor(icresboottestplot$Process,
                                         'Dispersal.Limitation'='DL',
                                         'Drift.and.Others'='DR',
                                         'Homogenizing.Dispersal'='HD',
                                         "Heterogeneous.Selection"="HeS",
                                         'Homogeneous.Selection'='HoS')




plot.icres.17.boot.bin24<-ggplot2::ggplot(data=icresboottestplot,
                                          aes(fill=Process,
                                              x=factor(Group,levels = c("172con",
                                                                        "172mon-6",
                                                                        '172mon1',
                                                                        '172mon8',
                                                                        '172mon15',
                                                                        '172mon22',
                                                                        '172mon50',
                                                                        '172mon78',
                                                                        '172mon106',
                                                                        '172mon134',
                                                                        '172mon372')),
                                              y=Observed))+
  geom_bar(position = "stack", stat="identity") +xlab("Days after injection") +
  theme_minimal() +ylab("Observed relative importance")

p.icres.17<-plot.icres.17.boot.bin24+scale_fill_manual(values=cbPalette[2:6])+
  scale_x_discrete(labels=c('Control','-6','1','8','15','22','50','78','106','134','372'))+
  theme(text = element_text( size = 20))

p.icres.17+ ggtitle("(b)")+theme(plot.title = element_text(hjust = 1.2, vjust=-3))
export::graph2ppt(file='icamp17')

ggsave(filename = "iCAMP17.jpg",
       plot=p.icres.17,
       path="~/R_files/ggsave_destination")
#####
#now icamp category
cate.17<-data.frame(matrix(NA, nrow=ntaxa(ps.EVO_2017_02_rar), ncol = 2))
cate.17[,1]<-taxa_names(ps.EVO_2017_02_rar)
cate.17[which(cate.17[,1] %in% taxa_names(ps.EVO_2017_02_rar_RA_stimactive)), 2]<-"stimulated"
cate.17[which(is.na(cate.17[,2])),2]<-"not_stimulated"
rownames(cate.17)<-cate.17[,1]
cate.17[,-1,drop=FALSE]->cate.17


clas.17<-data.frame(tax_table(ps.EVO_2017_02_rar))
icresbins.EVO_17<-icamp.bins(icres.17_bin24$detail, treat = NULL,
                             clas =clas.17 )

icrescate.EVO_17<-icamp.cate(icresbins.EVO_17,comm = icres.17_bin24$detail$comm,
                             cate = cate.17, treat = treat.icamp.17.boot)

icrescate.EVO_17$Ptx->iccate.17
iccate.17[,-c(1:2,9,15)]->iccate.17
iccate.17<-data.frame(pivot_longer(data=iccate.17, cols=2:11, names_to = c("Stimulation","Process"),
                                   names_sep = "[.]"))

plot.iccate.17<-ggplot2::ggplot(data=iccate.17, 
                                aes(fill=Process,
                                    x=factor(Group,levels = c("172con",
                                                              "172mon-6",
                                                              '172mon1',
                                                              '172mon8',
                                                              '172mon15',
                                                              '172mon22',
                                                              '172mon50',
                                                              '172mon78',
                                                              '172mon106',
                                                              '172mon134',
                                                              '172mon372')), 
                                    y=value)) +
  geom_bar(position = "stack", stat="identity")+xlab("Days after injection")+
  facet_grid(.~Stimulation, labeller= as_labeller(c('not_stimulated'="Not stimulated",'stimulated'='Stimulated')))+
  theme_minimal() +ylab("Observed relative importance")

p.icrescate.17<-plot.iccate.17+scale_fill_manual(values=cbPalette[2:6])+scale_fill_manual(values=cbPalette[2:6])+
  scale_x_discrete(labels=c('Control','-6','1','8','15','22','50','78','106','134','372'))+
  theme(text = element_text( size = 20))
p.icrescate.17+ ggtitle("(d)")+theme(plot.title = element_text(hjust = 1.2, vjust=-3),
                                     axis.text.x = element_text(size=14,angle=45,vjust = 1.5,  hjust=1.3))


#Figure_5:
get_legend(p.icres.09+theme(text = element_text( size = 12)))->icamp_legend

Fig_5<-gridExtra::grid.arrange(p.icres.09+
                                 coord_cartesian(clip='off')+
                                 theme(text = element_text( size = 12), legend.position = 'none',panel.grid.major.x = element_blank() ,
                                       plot.title = element_text(hjust = 1),axis.text.x = element_text(angle=90, hjust=1, vjust = 0.1))+
                                 ylab("Observed\nrelative importance")+
                                 ggtitle('a')+xlab(NULL),
                               
                               p.icres.17+
                                 coord_cartesian(clip='off')+
                                 theme(text = element_text( size = 12), legend.position = 'none',panel.grid.major.x = element_blank() ,
                                       plot.title = element_text(hjust = 1),axis.text.x = element_text(angle=90, hjust=1, vjust = 0.1))+
                                 ylab("Observed\nrelative importance")+
                                 ggtitle('b')+ylab(NULL)+xlab(NULL),
                               
                               icamp_legend,
                               
                               p.icrescate.09+
                                 coord_cartesian(clip='off')+
                                 
                                 theme(text = element_text( size = 12), legend.position = 'none',panel.grid.major.x = element_blank() ,
                                       plot.title = element_text(hjust = 1),axis.text.x = element_text(angle=90, hjust=1, vjust = 0.1))+
                                 ylab("Observed\nrelative importance")+
                                 ggtitle('c'),
                               
                               p.icrescate.17+
                                 coord_cartesian(clip='off')+
                                 
                                 theme(text = element_text( size = 12), legend.position = 'none',panel.grid.major.x = element_blank() ,
                                       plot.title = element_text(hjust = 1),axis.text.x = element_text(angle=90, hjust=1, vjust = 0.1, size=10))+
                                 ylab("Observed\nrelative importance")+
                                 ggtitle('d')+ylab(NULL),
                               #icamp_legend,
                               nrow=2,ncol=3,
                               widths=c(1,1,0.5)
)
ggsave('Figure_5.tiff', plot=Fig_5, device='tiff',width = 190, height =120, units='mm',dpi=1000 )

###############################################################################
#Running NST

treat.NST.09<-treat.icamp.09[,"phase", drop=FALSE]
save.wd.09.NST<-"~/R_files/EVO_09vs17_inputs_usearch/for_data_sharing/NST_outputs_09"
pd.big.09.NST=iCAMP::pdist.big(tree = tree.icamp.09, wd=save.wd.09.NST, nworker = 8)
pNSTres.09<-NST::pNST(comm=comm.icamp.09, 
                      tree=tree.icamp.09,
                      pd.wd=save.wd.09.NST,
                      group = treat.NST.09,
                      taxo.null.model = "PF",#null may be preferred
                      phylo.shuffle = TRUE,#w/ phylo shuffle
                      meta.group = NULL,
                      abundance.weighted = TRUE,
                      rand = 1000,
                      nworker = 8,
                      output.rand = TRUE)
pNSTbootres.09<-nst.boot(pNSTres.09,
                         group = treat.NST.09, 
                         out.detail = TRUE)

pNST09plotting<-as.data.frame(pivot_longer(as.data.frame(pNSTbootres.09$detail$NST.boot),cols=c(1:4), 
                                           names_to="Group",values_to="Obs" ))

p.pNST.09<-ggplot(data=pNST09plotting, mapping=(aes(x=Group,y=Obs))) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot()+
  theme_classic(base_size=15)+
  ylab('pNST')

nst.panova(pNSTres.09, group=treat.NST.09, rand=999, trace=TRUE)

tNSTres.09<-NST::tNST(comm=comm.icamp.09, 
                      group = treat.NST.09,
                      dist.method = "bray",
                      null.model = "PF",
                      abundance.weighted = TRUE,
                      rand = 1000,
                      nworker = 8,
                      output.rand = TRUE)

tNSTbootres.09<-nst.boot(tNSTres.09,
                         group = treat.NST.09, out.detail = TRUE)


tNST09plotting<-as.data.frame(pivot_longer(as.data.frame(tNSTbootres.09$detail$NST.boot),cols=c(1:4), 
                                           names_to="Group",values_to="Obs" ))

p.tNST.09<-ggplot(data=tNST09plotting, mapping=(aes(x=Group,y=Obs))) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot()+
  theme_classic(base_size=15)+
  ylab('tNST')

ggplot(data=tNST09plotting, aes(x=Group, y=Obs))+geom_violin(trim=F)+theme_classic()+
  xlab('tNST')+geom_boxplot(width=0.1)
ggplot(data=pNST09plotting, aes(x=Group, y=Obs))+geom_violin()+theme_classic()+
  xlab('pNST')





###############################################################################
#2017

treat.NST.17<-treat.icamp.17[,"phase", drop=FALSE]
save.wd.17.NST<-"~/R_files/EVO_09vs17_inputs_usearch/for_data_sharing/NST_outputs_17"
pd.big.17.NST=iCAMP::pdist.big(tree = tree.icamp.17, wd=save.wd.17.NST, nworker = 8)
pNSTres.17<-NST::pNST(comm=comm.icamp.17, 
                      tree=tree.icamp.17,
                      pd.wd=save.wd.17.NST,
                      group = treat.NST.17,
                      taxo.null.model = "PF",#null may be preferred
                      phylo.shuffle = TRUE,#w/ phylo shuffle
                      meta.group = NULL,
                      abundance.weighted = TRUE,
                      rand = 1000,
                      nworker = 8,
                      output.rand = TRUE)
pNSTbootres.17<-nst.boot(pNSTres.17,
                         group = treat.NST.17, 
                         out.detail = TRUE)
View(pNSTbootres.17$summary)
pNST17plotting<-as.data.frame(pivot_longer(as.data.frame(pNSTbootres.17$detail$NST.boot),cols=c(1:4), 
                                           names_to="Group",values_to="Obs" ))

p.pNST.17<-ggplot(data=pNST17plotting, mapping=(aes(x=Group,y=Obs))) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot()+
  theme_classic(base_size=15)+
  ylab('NST')

nst.panova(pNSTres.17, group=treat.NST.17, rand=999, trace=TRUE)

tNSTres.17<-NST::tNST(comm=comm.icamp.17, 
                      group = treat.NST.17,
                      dist.method = "bray",
                      null.model = "PF",
                      abundance.weighted = TRUE,
                      rand = 1000,
                      nworker = 8,
                      output.rand = TRUE)

tNSTbootres.17<-nst.boot(tNSTres.17,
                         group = treat.NST.17, out.detail = TRUE)
View(tNSTbootres.17$summary)

tNST17plotting<-as.data.frame(pivot_longer(as.data.frame(tNSTbootres.17$detail$NST.boot),cols=c(1:4), 
                                           names_to="Group",values_to="Obs" ))

p.tNST.17<-ggplot(data=tNST17plotting, mapping=(aes(x=Group,y=Obs))) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot()+
  theme_classic(base_size=15)+
  ylab('NST')

#####
#plots
p.tNST.09+ylim(0,1)+ylab("tNST")+
  scale_x_discrete(labels=c('Control','1','2','3'))+
  xlab('Injection phase')


p.tNST.17+ylim(0,1)+ylab('tNST')+
  scale_x_discrete(labels=c('Control','1','2','3'))+
  xlab('Injection phase')

p.pNST.09+ylim(0,1)+ylab('pNST')+
  scale_x_discrete(labels=c('Control','1','2','3'))+
  xlab('Injection phase')

p.pNST.17+ylim(0,1)+ylab('pNST')+
  scale_x_discrete(labels=c('Control','1','2','3'))+
  xlab('Injection phase')

#putting pNST and tNST results on the same plots, use geom_boxplot, colored by year

tNST09plotting$EVOYear<-c(rep("2009",nrow(tNST09plotting)))
tNST17plotting$EVOYear<-c(rep("2017",nrow(tNST17plotting)))
rbind(tNST09plotting, tNST17plotting)->tNSTplotting_combined
p.tNST<-ggplot(data=tNSTplotting_combined, mapping=(aes(x=Group,y=Obs))) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(aes(color=EVOYear))+
  theme_classic(base_size=15)+
  ylab('tNST')
p.tNST<-p.tNST+ylim(0,1)+
  scale_x_discrete(labels=c('Control','1','2','3'))+
  xlab('Injection phase')+
  theme(text = element_text(family = "sans", size = 15))+
  labs(color="Injection\nYear")

pNST09plotting$EVOYear<-c(rep("2009",nrow(pNST09plotting)))
pNST17plotting$EVOYear<-c(rep("2017",nrow(pNST17plotting)))
rbind(pNST09plotting, pNST17plotting)->pNSTplotting_combined
p.pNST<-ggplot(data=pNSTplotting_combined, mapping=(aes(x=Group,y=Obs))) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(aes(color=EVOYear))+
  theme_classic(base_size=15)+
  ylab('pNST')
p.pNST<-p.pNST+ylim(0,1)+
  scale_x_discrete(labels=c('Control','1','2','3'))+
  xlab('Injection phase')+
  theme(text = element_text(family = "sans", size = 15))+
  labs(color="Injection\nYear")

ggsave(filename = "tNST.jpg",
       plot=p.tNST,
       path="~/R_files/ggsave_destination")
ggsave(filename = "pNST.jpg",
       plot=p.pNST,
       path="~/R_files/ggsave_destination")



tNSTplotting_combined$NST<-rep('tNST', nrow(tNSTplotting_combined))
pNSTplotting_combined$NST<-rep('pNST', nrow(tNSTplotting_combined))
NSTplotting_combined<-rbind(tNSTplotting_combined,pNSTplotting_combined)

p.NST_09<-ggplot(data=subset(NSTplotting_combined,EVOYear==2009), aes(x=Group, y=Obs, fill=NST))+
  geom_violin(trim=T)+theme_classic()+
  xlab('Injection phase')+ylim(0,1.1)+ylab('Observed stochasticity ratio')+
  ggtitle('a')+theme(plot.title=element_text(hjust=1), legend.position='top')+
  scale_x_discrete(labels=c('Control' = 'Control', "X1" = "Phase 1",
                            "X2" = "Phase 2", 'X3'='Phase 3'))+
  labs(color='Normalized stochasticity ratio')

p.NST_17<-ggplot(data=subset(NSTplotting_combined,EVOYear==2017), aes(x=Group, y=Obs, fill=NST))+
  geom_violin(trim=F)+theme_classic()+
  xlab('Injection phase')+ylim(0,1.1)+ylab('Observed stochasticity ratio')+
  scale_x_discrete(labels=c('Control' = 'Control', "X1" = "Phase 1",
                            "X2" = "Phase 2", 'X3'='Phase 3'))+
  ggtitle('b')+theme(plot.title=element_text(hjust=1), legend.position='none')

Fig_S4<-gridExtra::grid.arrange(p.NST_09,
                                p.NST_17)
ggsave(plot = Fig_S4,
       filename = 'Figure_S4.tiff',
       device='tiff',
       dpi=1000,
       width=190,
       height=190,
       units='mm')

