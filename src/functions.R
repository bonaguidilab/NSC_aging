scalePerMillion = function (TPM){
  .d = apply(TPM, 2, function(X) return(X/sum(X))) * 1000000
  return (as.data.frame(.d))
}

pvolcanoPlotGG = function(log2FCs, pvals, genenames = NA, fcthresh = 2, alpha =.05, main = "Volcano Plot",highlight_cells=NULL,lft_nm = NULL,rt_nm=NULL,give_sig_names = FALSE){
  
  df <- data.frame(fc = log2FCs,
                   pv = -log10(pvals),
                   gene = genenames) %>% 
    mutate(significant = (pv >= (-log10(alpha)) & (abs(fc) >= fcthresh)))
  
  df$direction <- "NA"
  df$direction[df$significant & df$fc > 0] <- "UP"
  df$direction[df$significant & df$fc < 0] <- "DOWN"
  
  
  df_sig <- df %>% 
    filter(significant)
  
  df_notsig <- df %>%
    filter(!significant)
  
  
  p <- ggplot()+
    geom_hline(yintercept = -log10(alpha),linetype = "dashed",size = .25)+
    geom_vline(xintercept = fcthresh,linetype = "dashed",size = .25)+
    geom_vline(xintercept = -fcthresh,linetype = "dashed",size = .25)+
    geom_point(mapping = aes(fc,pv),
               data = df_notsig,
               size = .25)+
    geom_point(mapping = aes(fc,pv,color = direction),
               data = df_sig,
               size = 2,
               inherit.aes = FALSE)+
    scale_x_continuous(limits = c(-max(abs(df$fc)),max(abs(df$fc))))+
    labs(title = main,
         x = expression("log"[2]*"(Fold Change)"),
         y = expression("-log"[10]*"(adj-p)"),
         color = "Direction")+
    scale_color_brewer(palette = "Set1",direction = -1)+
    theme_bw()+
    theme(plot.title = element_text(hjust = .5,size = 15))
  
  if (give_sig_names){
    p <- p %>% ggrepel::geom_label_repel(mapping = aes(fc,pv,color = direction,color = direction,label = gene),
                                data = df_sig,
                                fill = scales::alpha("grey",.25),
                                inherit.aes = FALSE)
  }
  
  
  
  candidates = list(upGenes = as.character(df_sig$gene)[df_sig$direction=="UP"], dnGenes = as.character(df_sig$gene)[df_sig$direction=="DOWN"],plot = p) 
  
  return(candidates)
}


expandEdgeCols = function(o.mat,gene_cols){
  edge_cols = c()
  for (i in 1:nrow(o.mat)){
    count = sum(as.numeric(o.mat[i,]))
    edge_cols = c(edge_cols,rep(gene_cols[i],count))
  }
  return(edge_cols)
}


getGeneOccupancyGroup = function(o.mat, seed = 5, use_tSNE = FALSE,col_options = NULL){
  
  gene_vect = rep(0,nrow(o.mat))
  names(gene_vect) = rownames(o.mat)
  
  type = 1
  for(i in 1:nrow(o.mat)){
    if (gene_vect[i]==0){
      tmplt = as.numeric(o.mat[i,])
      gene_vect[i] = type
      
      for (h in which(gene_vect==0)){
        query = as.numeric(o.mat[h,])
        if (all(tmplt == query)){
          gene_vect[h] = type
        }
        
      }
      type = type +1
    }
  }
  
  
  set.seed(seed)
  
  n = length(unique(gene_vect))
  
  #if (n > 40) use_tSNE = TRUE
  if(is.null(col_options)){
    color = distinctColorPalette(k = n, altCol = FALSE, runTsne = use_tSNE)
    
  }else{
    color = col_options
  }
  gene_vect[names(gene_vect)] = color[as.numeric(gene_vect)]
  return(gene_vect)
}



get_sig_terms_return_gd = function(up,dn, gene_universe, p_adj = TRUE){
  up_raster = factor(as.integer(gene_universe %in% up))
  names(up_raster) = gene_universe
  
  dn_raster = factor(as.integer(gene_universe %in% dn))
  names(dn_raster) = gene_universe
  
  ontgs = c("MF","BP","CC")
  up = c()
  dn = c()
  gd_lst_up = list()
  gd_lst_dn = list()
  for (ont in ontgs){
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    
    up_GD = new("topGOdata", ontology = ont, allGenes = up_raster, geneSel = function(p) p < 
                  0.9, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", 
                ID = "symbol")
    dn_GD = new("topGOdata", ontology = ont, allGenes = dn_raster, geneSel = function(p) p < 
                  0.9, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", 
                ID = "symbol")
    
    gd_lst_up[[ont]] = up_GD
    gd_lst_dn[[ont]] = dn_GD
    up_rf = getSigGroups(up_GD, test.stat)
    dn_rf = getSigGroups(dn_GD, test.stat)
    
    up = c(up,score(up_rf))
    dn = c(dn,score(dn_rf))
  }
  
  up_all = c()
  up_in = c()
  dn_all = c()
  dn_in = c()
  for (ont in ontgs){
    
    ug = genesInTerm(gd_lst_up[[ont]],names(up)[which(Ontology(names(up))==ont)])
    dg = genesInTerm(gd_lst_dn[[ont]],names(dn)[which(Ontology(names(dn))==ont)])
    
    up_all = c(up_all,
               unlist(lapply(ug,length)))
    up_in = c(up_in,
              unlist(lapply(ug,function(v){
                return(length(intersect(v,names(up_raster[up_raster==1]))))
              })))
    dn_all = c(dn_all,
               unlist(lapply(dg,length)))
    dn_in = c(dn_in,
              unlist(lapply(dg,function(v){
                return(length(intersect(v,names(dn_raster[dn_raster==1]))))
              })))
    
  }
  
  up_all = up_all[names(up)]
  up_in  = up_in[names(up)]
  dn_all = dn_all[names(dn)]
  dn_in  = dn_in[names(dn)]
  
  
  if (p_adj){
    up = p.adjust(up,method = "BH")
    dn = p.adjust(dn,method = "BH")
  }
  
  
  return(list(up = list(up = up,gd = gd_lst_up,up_all = up_all, up_in = up_in ),
              dn = list(dn = dn,gd = gd_lst_dn,dn_all = dn_all, dn_in = dn_in)))
}



GOstringNetwork = function(gd,score,ontology = "BP",thresh = .05,
                           main_title = "",
                           color_by_connectivity = FALSE,
                           take_top = NULL,
                           min_connections = 1,
                           mode = c("kamadakawai","circle","fruchtermanreingold","target"),
                           seed = 5,
                           alpha = 1,
                           constrained_to = NULL,
                           height = 100,
                           width = 100,
                           padding = .25,
                           gene_text = 5,
                           term_text = 8.9,
                           gene_size = 1,
                           term_size = 5,
                           paint_edges_as_nodes = FALSE,
                           edge_cols = "grey50",
                           actor_col = "grey",
                           event_col = "gold",
                           background = "black",
                           force = 1){
  
  if (length(mode)==4){
    mode = mode[1]
  }
  
  
  
  sig_terms = sort(score)
  sig_terms = names(sig_terms)[sig_terms<thresh]
  sig_terms = sig_terms[Ontology(sig_terms)==ontology]
  
  if ( (!is.null(constrained_to)) & (!is.null(take_top)) ){
    stop("Arguments for both 'constrained_to' and 'take_top' were non NULL")
  }
  
  if (!is.null(constrained_to)){
    sig_terms = sig_terms[sig_terms %in% constrained_to]
    
    not_in = constrained_to[!constrained_to %in% sig_terms]
    if (length(not_in) >0){
      warning(paste0(length(not_in)," terms in 'constrained_to' not present in score \n ",not_in))
    }
    
  }
  
  if (length(sig_terms) <1){
    stop("constraints on 'score' leave no remaining terms")
  }
  
  
  if (!is.null(take_top)){
    sig_terms = sig_terms[1:take_top]
  }
  
  genes_in_sig_terms = genesInTerm(gd, sig_terms)
  
  census_genes = sigGenes(gd)
  census_genes = census_genes[census_genes %in% unique(unlist(genes_in_sig_terms))]
  
  o.mat = matrix(0,nrow = length(census_genes),ncol = length(genes_in_sig_terms))
  rownames(o.mat) = census_genes
  colnames(o.mat) = names(genes_in_sig_terms)
  o.mat = as.data.frame(o.mat)
  
  for (term in names(genes_in_sig_terms) ){
    o.mat[,term] = as.numeric(census_genes %in% genes_in_sig_terms[[term]])
  }
  
  keep_rows_bool = apply(o.mat,1,function(x){return(sum(x)>=min_connections)})
  discarded_genes = rownames(o.mat)[which(!keep_rows_bool)]
  
  o.mat = o.mat[which(keep_rows_bool),]
  
  keep_cols_bool = apply(o.mat,2,function(x){return(sum(x)>0)})
  discarded_terms = colnames(o.mat)[which(!keep_cols_bool)]
  
  o.mat = o.mat[,which(keep_cols_bool)]
  
  if (length(discarded_genes >0)){
    cat(paste("Keeping",length(which(keep_rows_bool)),"genes and discarding",length(discarded_genes),"genes that have fewer than",as.character(min_connections),"connections \n"))
    if(length(discarded_terms)>0){
      cat(paste("Discarding",length(discarded_terms),"because they have 0 connections \n"))
    }
  }
  
  
  
  
  f.mat = o.mat
  colnames(f.mat) = Term(colnames(o.mat))
  
  bip = network(f.mat,
                matrix.type = "bipartite",
                ignore.eval = FALSE,
                names.eval = "weights")
  
  gene_cols = rep(actor_col,nrow(f.mat))
  if(color_by_connectivity) gene_cols = getGeneOccupancyGroup(f.mat, seed = seed)
  if(paint_edges_as_nodes) edge_cols = expandEdgeCols(f.mat, gene_cols)
  gene_cols = alpha(gene_cols,alpha)
  
  set.seed(1)
  p <- ggnet2(bip,shape = "mode",color = "mode", 
         node.color = c(gene_cols,
                        rep(event_col,ncol(f.mat))), 
         edge.color =  edge_cols,
         size = c(rep(gene_size,nrow(f.mat)),
                  rep(term_size,ncol(f.mat))),
         label = FALSE,mode = mode) +
    ggtitle(main_title)+
    geom_text_repel(aes(label = c(rownames(f.mat),
                                  colnames(f.mat) )),
                    cex = c(rep(gene_text,nrow(f.mat)),
                            rep(term_text,ncol(f.mat))),
                    min.segment.length = unit(.2, "lines"))+
    theme(legend.position="none",plot.title = element_text(hjust = .5))

  return(p)
  
}


#################
# Velocity stuff
#################
grp_avg_origin = function(all_coords,group_lst = list()){
  grps = names(group_lst)
  avg_mat = matrix(0,nrow = length(grps),ncol = ncol(all_coords))
  rownames(avg_mat) = grps
  colnames(avg_mat) = c("x","y")
  
  for (grp in grps){
    cls = group_lst[[grp]]
    
    x_avg = mean(all_coords[cls,1])
    y_avg = mean(all_coords[cls,2])
    
    avg_mat[grp,] = c(x_avg,y_avg)
  }
  return(avg_mat)
}

get_transposed_vectors = function(all_coords,group_lst = list(),scale = 1){
  grps = names(group_lst)
  
  grp_origin = grp_avg_origin(all_coords[,1:2],group_lst)
  
  oc_coords = all_coords - cbind(all_coords[,1:2],all_coords[,1:2])
  oc_coords = t(t(oc_coords)* c(1,1,scale,scale))
  
  for (grp in grps){
    cls = group_lst[[grp]]
    ogn = grp_origin[grp,]
    
    oc_coords[cls,1] = rep(ogn[1],length(cls))
    oc_coords[cls,2] = rep(ogn[2],length(cls))
    
    oc_coords[cls,3] = oc_coords[cls,3] + ogn[1]
    oc_coords[cls,4] = oc_coords[cls,4] + ogn[2]
  }
  
  
  return(oc_coords)
  
}

get_average_group_vector = function(all_coords,group_lst = list(),comparison_lst = list(), scale = 1){
  grp_origin = grp_avg_origin(all_coords[,1:2],group_lst)
  transposed_coords = get_transposed_vectors(all_coords,group_lst,scale = scale)
  
  grps = names(group_lst)
  cmps = names(comparison_lst)
  
  comp_mats = list()
  for (cmp in cmps){
    tmp_cmp_mat = matrix(0,nrow = length(grps),ncol = ncol(all_coords))
    rownames(tmp_cmp_mat) = grps
    colnames(tmp_cmp_mat) = colnames(all_coords)
    
    for (grp in grps){
      cls = intersect(group_lst[[grp]],comparison_lst[[cmp]])
      
      tmp_cmp_mat[grp,] = apply(transposed_coords[cls,],2,mean)
      
    }
    comp_mats[[cmp]] = tmp_cmp_mat
  }
  
  
  return(comp_mats)
}

