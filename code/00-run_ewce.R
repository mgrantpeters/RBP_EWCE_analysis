library("EWCE")
library("biomaRt")
library("ggplot2")
library("grid")
set.seed(1234)
setwd("/Users/melis/Documents/SoniaGR_RBP_EWCE/code/")

load("C:/Users/melis/Downloads/ctd_aibsMultipleCrtxSmrtSeq.rda")  #Load CTD data
n = 10000

data = read.csv("../data/RBPs_subgroups.csv")

list = data[data$Spliceosome>0,1]

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 1,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length
thisResult1 = data.frame(full_results$results)
thisResult1$MajorCellType = unlist(lapply(strsplit(as.character(thisResult1$CellType),'_'), `[[`, 1))
thisResult1$testList = 'Spliceosome'
FinalResult1 = thisResult1


plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "FDR",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="../plots/EWCE_RBPs_spliceosome_annot1.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "../processed-data/EWCE_RBPs_spliceosome_annot1.csv", row.names=TRUE, quote=FALSE) 

rm(full_results)
rm(plot_list)

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 2,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length

thisResult = data.frame(full_results$results)
thisResult$MajorCellType = unlist(lapply(strsplit(as.character(thisResult$CellType),'_'), `[[`, 1))
thisResult$testList = 'Spliceosome'
FinalResult = thisResult

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "FDR",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="../plots/EWCE_RBPs_spliceosome_annot2.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 20,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "../processed-data/EWCE_RBPs_spliceosome_annot2.csv", row.names=TRUE, quote=FALSE) 
rm(full_results)
rm(plot_list)

#################------------------- Splicing regulation
#Load target gene list

list = data[data$Splicing.regulation>0,1]

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 1,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length

thisResult1 = data.frame(full_results$results)
thisResult1$MajorCellType = unlist(lapply(strsplit(as.character(thisResult1$CellType),'_'), `[[`, 1))
thisResult1$testList = 'Splicing regulation'
FinalResult1 = rbind(FinalResult1,thisResult1)

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "FDR",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="../plots/EWCE_RBPs_splicing-regulation_annot1.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "../processed-data/EWCE_RBPs_splicing-regulation_annot1.csv", row.names=TRUE, quote=FALSE) 

rm(full_results)
rm(plot_list)

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 2,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length
thisResult = data.frame(full_results$results)
thisResult$MajorCellType = unlist(lapply(strsplit(as.character(thisResult$CellType),'_'), `[[`, 1))
thisResult$testList = 'Splicing regulation'
FinalResult = rbind(FinalResult,thisResult)

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "FDR",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="../plots/EWCE_RBPs_splicing-regulation_annot2.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 20,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "../processed-data/EWCE_RBPs_splicing-regulation_annot2.csv", row.names=TRUE, quote=FALSE) 
rm(full_results)
rm(plot_list)

############-------------------------------- FINAL SAVE and PLOT

write.csv((.packages()), "../processed-data/loaded_R_packages.R")


write.csv(FinalResult1, "../processed-data/final_results_table_annot1.csv")
sig1 = FinalResult1[FinalResult1$q>0.05,]

   
FinalResult1 %>% 
  mutate(CellType = factor(CellType,c("Astrocyte", "Microglia","OPC", "Oligodendrocyte","Vascular_cells","GABAergic_LAMP5","GABAergic_PAX6","GABAergic_PVALB","GABAergic_SST","GABAergic_VIP","Glutamatergic_IT","Glutamatergic_L4_IT","Glutamatergic_L5_ET","Glutamatergic_L5_6_IT_Car3","Glutamatergic_L5_6_NP", "Glutamatergic_L6_CT","Glutamatergic_L6b"))) %>%        
  ggplot() +
  geom_point(aes(x = -log10(q), y = CellType,color = MajorCellType, size = fold_change, alpha = q<0.05))+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(rows = vars(testList))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(q)") +
  ylab("Cell type") +
  coord_flip() + # this is to make  the graph landscape
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  guides(fill = guide_legend(title = "Major Cell Type")) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=90, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.title = element_text(colour = "black", size = "10"),
        legend.position = "top",
        ## This is to plot the two legends in two rows
        legend.box="vertical") 
ggsave(
  filename="../plots/Overview_results_annot1_test.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 8,
  height = 7,
  units = c("in"),
  dpi = 300
)

sig = FinalResult[FinalResult$p<0.05,]

write.csv(FinalResult, "../processed-data/final_results_table_annot2.csv")

ggplot(FinalResult) +
  geom_point(aes(x = -log10(q), y = CellType, color = MajorCellType, size = fold_change, alpha = q<0.05))+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(rows = vars(testList))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(q)") +
  ylab("Cell type") +
  coord_flip() + # this is to make  the graph landscape
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  guides(fill = guide_legend(title = "Major Cell Type")) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=90, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.title = element_text(colour = "black", size = "10"),
        legend.position = "top",
        ## This is to plot the two legends in two rows
        legend.box="vertical") 

ggsave(
  filename="../plots/Overview_results_ann2.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 22,
  height = 15,
  units = c("in"),
  dpi = 300
)