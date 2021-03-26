source("src/functions.R") #collection of custom functions for processing and visualization
source("src/heatmap3.R")

library("tidyverse")
library("GEOquery")
library("GO.db")
library("org.Mm.eg.db")
library("topGO")
library("edgeR")
library("gridExtra")
library("RColorBrewer")
library("ggnet")
library("network")
library("randomcoloR")
library("ggrepel")


#load data
ap.cts <-  read.csv("data/expression_mats_counts.csv",row.names = 1,header = TRUE)
ap.cpm <- scalePerMillion(ap.cts)
ap.tpm <- read.csv("data/expression_mats_TPM.csv",row.names = 1,header = TRUE) #assign 'new = TRUE' to generate new expression table
ap.grps <- read.csv("data/expression_mats_groups.csv",row.names = 1,header = TRUE)

ap.pca <- readRDS(file = "data/magdif/max_pca_num.rds"); ap.pca = ap.pca[colnames(ap.tpm),]
ap.cols = readRDS(file = "data/magdif/colors.rds"); ap.cols = ap.cols[colnames(ap.tpm)]


#cell state groups 
qnsc_cls <- as.character(ap.grps$names[ap.grps$groups %in% c("qNSC") ])
ansc_cls <- as.character(ap.grps$names[ap.grps$groups %in% c("aNSC") ])
ipc_cls <- as.character(ap.grps$names[ap.grps$groups %in% c("IPC") ])


###############################################
###  FIGURE 3C & SUPPLEMENTAL FIGURE S4A-C  ###
###############################################
type_gns = list()
type_gns[["NSC"]] = c("Ntsr2","Nes","Vim","Fabp7")
type_gns[["Quiescent"]] <- c("Aldoc","Aqp4","Id3","Hes1")
type_gns[["Proliferating"]] = c("Mki67","Pcna","Mcm6")
type_gns[["IPC"]] = c("Eomes","Stmn1","Sox11")


par(mar = c(3,3,3,3),mfrow = c(2,2))
cex = .5

length = 100
marker = colorRampPalette(brewer.pal(11,"Spectral"))(length)[length:1]
for (tp_gn in names(type_gns)){
  
  gns = type_gns[[tp_gn]]
  
  #exp = apply(ap.tpm[gns,],2,gm_mean)
  exp = log2(apply(ap.tpm[gns,],2,gm_mean)+1)
  exp_col = exp/max(exp)
  exp_col = exp_col + .01
  
  clrs = marker[as.numeric(round(exp_col*length,0))]
  plot(ap.pca[,1:2],cex = cex,xlab = "PC1 1",ylab = "PC 2",main = tp_gn)
  points(ap.pca[,1:2],cex = cex-(.1*cex),col = alpha(clrs, 1),pch = 19)
  
  legend.col(col = marker, lev = as.numeric(exp))
}



###############################################
###    DERIVATION OF DEX GENES WITH edgeR  ###
###############################################

#edgeR derivation of differentially expressed genes between young qNSCs and old qNSCs. 
x = ap.cts[,qnsc_cls] # assign count matrix to variable 
group = rep(0,ncol(x)) #initialize group vector, elements correspond to columns across the count matrix
group[grep("young",colnames(x))] = 1 # group 1 is assigned to indices that correspond to YOUNG cells in count matrix
group[grep("old",colnames(x))] = 2 #group 2 is assigned to indices that correspond to OLD cells in count matrix 
group = factor(group) # give group vector factors (this is necessary for edgeR)

#feature selection: this loop simply finds the features (genes) which have a mean count of greater than 5 in young or old cells
kg = apply(x,1,function(v){
  keep_it = FALSE
  if (mean(v[group ==1])>5|mean(v[group ==2])>5 ){
    keep_it = TRUE
  }
  return(keep_it)
})

#include only the kepg genes 
x = x[kg,]

y = DGEList(counts=x,group=group) #creates data structure for donstream functions
y = calcNormFactors(y) #calculates normalization factors to scale the raw library sizes 
design = model.matrix(~group) #creates a design matrix by expanding factors to a set of dummy variables
y = estimateDisp(y,design) #Maximizes the negative bionomial likelihood to give the estimate of common, trended and tagwise dispersions
#Perform quasi-likelihood F-tests:
fit = glmQLFit(y,design) # fit to a quasi-likelihood negative binomial generalized log-linerar model to count data. 
qlf = glmQLFTest(fit,coef=2) #Conduct genewise test. verbose information about the piepeline is all assigned to this variable
qlf$table$FDR = p.adjust(qlf$table$PValue,method = "BH") #insert FDR to qlf data frame 


#Biological Coefficient of Variation
###############################
### SUPPLEMENTAL FIGURE S4D ###
###############################
plotBCV(y)

# volcano plot which sets FDR and FC thresholds
old_vcp_nsc_edg = pvolcanoPlotGG(qlf$table$logFC,p.adjust(qlf$table$PValue,method = "BH"),
                                 rownames(qlf$table),
                                 alpha = .05,
                                 main = "Volcano - Old/Young [edgeR]")

#visualize volcano plot
###############################
### SUPPLEMENTAL FIGURE S4E ###
###############################
plot(old_vcp_nsc_edg$plot)

###################
###  FIGURE 3G  ###
###################
exp = qlf$fitted.values[,qnsc_cls]

young_cls = grep("young",colnames(exp))
old_cls = grep("old",colnames(exp))
exp = log2(exp[,c(young_cls,old_cls)]+1)
column_cols = t(t(c(rep("#0060C3",length(young_cls)),
                    rep("#ff4d4d",length(old_cls)))))
row_cols = t(c(rep("#0060C3",length(old_vcp_nsc_edg$dnGenes)),
                    rep("#ff4d4d",length(old_vcp_nsc_edg$upGenes))))
mat = as.matrix(exp[c(old_vcp_nsc_edg$dnGenes,old_vcp_nsc_edg$upGenes),])


v = heatmap.3(mat,
              col = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                      "RdYlBu")))(100),
              #add_grid = FALSE,
              ColSideColors = column_cols,
              #RowSideColors = row_cols,
              KeyValueName = expression('log'[2]*"(TPM+1)"),
              Rowv = NA,
              Colv = NA,
              cexCol = .0001,
              cexRow = .0001)
legend("topright", legend=c("2mo", "4.5mo"),
       col=c("#0060C3","#ff4d4d"), pch = 19, cex=0.8)



#####################
###  FIGURE 3J-L  ###
#####################

#these are the genes detected in the enriched go terms 
exp_lst = list("Generation of neurons" = c("Igf1", "Sox11","Trim11","Wnt5b" ),
               "Positive regulation of neurogenesis" = c("Socs2", "Wnt3", "Epha4", "Zfp365", "Wnt5a"), #expected
               "Gliogenesis" = c("Ezh2", "Disc1", "Mag", "Plp1", "Aspa", "P2ry12", "Tnfrsf21" ),
               "Cell cycle" = c("Heca","Wee1","Mcm6","Nup43", "Nsl1") #expected
               )

sig_lst = list("Negative regulation of signaling" = c("Lef1", "Sesn2", "Nfkbil1", "Abl1", "Abl2", "Crh"),
               "Semaphorin signaling" =c("Plxna4", "Plxnb3", "Nrp2", "Farp2"), #signaling
               "Ras signaling" = c("Rassf1","Icmt","Arhgap32","Cdh13","Arhgef11"), 
               "Rho signaling" = c("Spata13", "Myo9b", "Tiam2", "Mcf2l", "Scai") #signaling
               )

age_lst = list("Histone demethylation" = c("Kdm1b", "Kdm2a", "Kdm4a", "Kdm5d", "Uty"), #aging
               "Transcription" = c("Mycn", "Creb1", "Hif1an", "Sertad2", "Epcam","Sepsecs","Naf1","Sox30"),
               "NIK|NF kappaB signaling" =c("Cd14", "Malt1", "Eif2ak2", "Tirap", "Myd88"), 
               "DNA recombination" = c("Atad5", "Hus1", "Pms2", "Lig3", "Mcm9","Tep1"),
               "DNA repair" = c("Msh3", "Msh6", "Rfc3", "Rfc5", "Rad18","Pold1"), #TODO: left off ehre 
               "Double strand break repair" = c("Chek1", "Blm", "Spidr", "Fancb") # aging
               )

lsts = c("exp_lst","sig_lst","age_lst")

exp = qlf$fitted.values[,qnsc_cls]
young_cls = grep("young",colnames(exp))
old_cls = grep("old",colnames(exp))

exp = log2(exp[,c(young_cls,old_cls)]+1)
mat = as.matrix(exp[c(old_vcp_nsc_edg$dnGenes,old_vcp_nsc_edg$upGenes),])
column_cols = t(t(c(rep("#0060C3",length(young_cls)),
                    rep("#ff4d4d",length(old_cls)))))
 
for (k in 1:length(lsts)){
  lst = lsts[k]
  tmp_lst = get(lst)
  
  tmp_lst = lapply(tmp_lst,function(v){
    return(sort(v))
  })
  
  gns_vct = as.vector(as.character(unlist(tmp_lst)))

  mat2 = mat[gns_vct,]
  
  set.seed(878+k)
  clrs = randomcoloR::distinctColorPalette(k = length(tmp_lst),runTsne = TRUE)
  row_cols = c()
  for (i in 1:length(tmp_lst)){
    row_cols = c(row_cols,rep(clrs[i],length(tmp_lst[[i]])))
  }
    
  row_cols = t(row_cols)
  print(row_cols)
  heatmap.3(mat2,
            col = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                    "RdYlBu")))(100),
            ColSideColors = column_cols,
            RowSideColors = row_cols,
            KeyValueName = "Log2(edgeR+1)",
            Rowv = NA,
            Colv = NA,
            cexRow = 1,
            cexCol = .001)
  
  legend("topright", legend=names(tmp_lst),
         col=clrs, pch = 19, cex=.75,ncol = 3)

}



#get qNSC associated GO terms
qNSC_trms = get_sig_terms_return_gd(old_vcp_nsc_edg$upGenes, old_vcp_nsc_edg$dnGenes,rownames(ap.cts),p_adj = FALSE)

##################
###  FIGURE 3H ###
##################
sig_terms_up <- qNSC_trms$up$up < .05
up_df <- tibble(GO_name = Term(names(qNSC_trms$up$up[sig_terms_up])),
                GO_id = names(qNSC_trms$up$up[sig_terms_up]),
                ontology = Ontology(names(qNSC_trms$up$up[sig_terms_up])),
                p_value = qNSC_trms$up$up[sig_terms_up],
                n_sig_genes = qNSC_trms$up$up_in[sig_terms_up],
                n_total_genes = qNSC_trms$up$up_all[sig_terms_up])

##################
###  FIGURE 3I ###
##################
sig_terms_dn <- qNSC_trms$dn$dn < .05
dn_df <- tibble(GO_name = Term(names(qNSC_trms$dn$dn[sig_terms_dn])),
                GO_id = names(qNSC_trms$dn$dn[sig_terms_dn]),
                ontology = Ontology(names(qNSC_trms$dn$dn[sig_terms_dn])),
                p_value = qNSC_trms$dn$dn[sig_terms_dn],
                n_sig_genes = qNSC_trms$dn$dn_in[sig_terms_dn],
                n_total_genes = qNSC_trms$dn$dn_all[sig_terms_dn])



###########################################
### Figure 4A & SUPPLEMENTAL FIGURE S4G ###
###########################################

up_trm_lst = list(post_translational_modifications = c("GO:0042158","GO:0006464"),
                  dna_repair = c("GO:0045739","GO:2000779"),
                  cell_signaling = c("GO:0071526","GO:1901224"),
                  cell_morphogenesis = c("GO:0050773","GO:0008366"),
                  neurogenesis = c("GO:0022008"),
                  giogenesis = c("GO:0042063"))

p_up = GOstringNetwork(qNSC_trms$up$gd$BP,qNSC_trms$up$up,ontology = "BP",thresh = .05,
                       main_title = paste0("UP Network min ",1," connections"),
                       #mode = "fruchtermanreingold",
                       constrained_to = as.vector(unlist(up_trm_lst)),
                       color_by_connectivity = TRUE,
                       paint_edges_as_nodes = TRUE,
                       min_connections = 1,
                       seed = 243453,
                       alpha = .7,
                       height = 25,
                       width = 25,
                       gene_text = 5,
                       term_text = 8.9,
                       gene_size = .25,
                       term_size = 1.25,
                       force = 2,
                       actor_col = alpha("#ff4d4d",.4))
plot(p_up) #try using ggsave() to make it less cluttered


###############################
### SUPPLEMENTAL FIGURE S4G ###
###############################
down_trm_lst = list(epigenetics = c("GO:0010608"), #"GO:0006464"
                    transcription = c("GO:0010468"),
                    dna_repair = c("GO:0006281"), #"GO:0006310"
                    cell_cycle = c("GO:0007049"),
                    cell_signaling = c("GO:0023057","GO:0007265"),
                    differentiation = c("GO:0022008"))

p_down = GOstringNetwork(qNSC_trms$dn$gd$BP,qNSC_trms$dn$dn,ontology = "BP",thresh = .05,
                               main_title = paste0("DOWN Network min ",1," connections"),
                               mode = "fruchtermanreingold",
                               constrained_to = as.vector(unlist(down_trm_lst)),
                               color_by_connectivity = TRUE,
                               paint_edges_as_nodes = TRUE,
                               min_connections = 1,
                               seed = 243453,
                               alpha = .7,
                               height = 25,
                               width = 25,
                               gene_text = 5,
                               term_text = 8.9,
                               gene_size = 1/2,
                               term_size = 5/2,
                               force = 2,
                               actor_col = alpha("#ff4d4d",.4))

plot(p_down)

##################
###  FIGURE 4B ###
##################
qnsc_exp = log2(ap.tpm[,qnsc_cls]+1)

#make batplots 
qnsc_mat_tmp = matrix(as.matrix(qnsc_exp),
                 ncol = ncol(qnsc_exp),
                 nrow = nrow(qnsc_exp),
                 dimnames = list(gene_name = rownames(qnsc_exp),
                                 cells = colnames(qnsc_exp)))
df_qnsc = as.data.frame(as.table(qnsc_mat_tmp))
df_qnsc$age = rep("2_mo",nrow(df_qnsc))
df_qnsc$age[grep("old",df_qnsc$cells)] = "4.5_mo"


df_qnsc_select = df_qnsc[df_qnsc$gene_name %in% c("Abl1","Abl2"),]

p_abl <- ggplot(df_qnsc_select,aes(x = gene_name,y = Freq, fill = age)) + 
  geom_violin(position=position_dodge(1),
              scale = "width")+
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width = .15,
                                                      dodge.width = 1),aes(colour = age))+
  scale_fill_manual(values=scales::alpha(c("#f8b195","#355c7d"),.5))+
  scale_color_manual(values =c("#f8b195","#355c7d")  )+
  ggtitle("Alb1 & Alb2 Expression in qNSCs")+
  ylab("log2(TPM+1)")+
  xlab("Gene")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plot(p_abl)

