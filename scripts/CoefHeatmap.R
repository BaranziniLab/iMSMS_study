# @diff.species a list of differential species in multiple comparisons from {linear_regression} function.
# @selected a vector of selected feature to plot, default NULL
# @name specify the name of column with p or fdr value
# @p.value cutoff of p-value to select the feature to plot. 
# @fdr cutoff of fdr to select the feature to plot.
# @coef boolean values determining if features should be selected by coefficient
# @colnames column names shown in heatmap
# @angle_col angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315)
# @cluster_rows,cluster_cols boolean values determining if rows/columns should be clustered
# @breaks boolean values determining if breaks of legend will be specified
# @width/height width/height of heatmap
# @out.file  name of output file


require(pheatmap)
require(RColorBrewer)

CoefHeatmap = function(diff.species, selected =NULL, name  =NULL,p.value=0.05, fdr = 0.05,coef =FALSE,
                       colnames = c("MS", "RRMS","PMS"),angle_col= 0,cluster_rows =T,cluster_cols=F, breaks =FALSE, width = 8, height =10,border.col="grey",
                       out.file = "results/WOL_shogun/Linear_regression/New202104/TreatedMS_RRMS_PMS_HHC_heatmap_fdr0.05.pdf"){
  
# get the p or fdr column index 
    if(!is.null(fdr)){
      if(!is.null(name)){
        index = paste("fdr",name,sep="_")
      }else{
        index = 17
      }
      if(is.null(selected )){
      sps = lapply(diff.species, function(x){
        y  = x[x[,index] <= fdr, ]
        y$taxonomy
      })
      }
    }else{
      if(!is.null(name)){
        index = paste("Pr",name,sep="_")
      }else{
        index = 5
      }
      if(is.null(selected )){
      sps = lapply(diff.species , function(x){
        y  = x[x[,index] <= p.value, ]
        if(coef){
          quantiles = quantile(x[,1], seq(0,1,0.05))
          min = quantiles["5%"]
          max = quantiles["95%"]
          y= y[y[,1] <= min | y[,1] >= max, ]
        }
        y$taxonomy
      })
      }
    }
    
  if(is.null(selected )){
  sps =unique(unlist(sps))
  }else{
    sps = selected
  }
  
    data = lapply(diff.species ,function(x){
      if(!is.null(name)){
        y = x[x$taxonomy %in% sps, c("taxonomy", paste("Coef", name, sep="_"),index)]
      }else{
        
        y = x[x$taxonomy %in% sps, c(which(colnames(x) == "taxonomy"),1,index)]
      }
      y
    })
  
  dataall = Reduce(function(x,y) merge(x,y, by = "taxonomy",all =T),data )
  dataall= dataall[match(sps, dataall$taxonomy),]
  #dataall[is.na(dataall)] = 0
  rownames(dataall) = dataall$taxonomy
  dataall = dataall[,-1]
  coefs = dataall[,grepl("Coef", colnames(dataall))]
  if(!is.null(fdr)){
    pvls = dataall[,grepl("fdr", colnames(dataall))]
  }else{
    pvls = dataall[,grepl("Pr", colnames(dataall))]
  }
  colnames(coefs) =colnames(pvls) = colnames
  coefs[is.na(coefs)] = 0
  pvls[is.na(pvls)] = 1
  pvls = ifelse(pvls < 0.001, "***", ifelse(pvls < 0.01, "**",ifelse(pvls < 0.05, "*", ifelse(pvls ==1,"NA", ""))))
  #my.breaks <- c(seq(min(coefs), max(coefs), by=0.01))
  
  coefs = coefs[!grepl("UNMAPPED|UNINTEGRATED|Unclassified|Unknown", rownames(coefs)),]
  pvls = pvls[!grepl("UNMAPPED|UNINTEGRATED|Unclassified|Unknown", rownames(pvls)),]
  
  coef.min = min(coefs)
  coef.max = max(coefs)
  if(abs(coef.min) > coef.max){
    max = abs(coef.min)
  }else{
    max = coef.max
  }
if(breaks){
  breaksList1 = seq(coef.min, 0,by = 0.01)
  breaksList2 = seq(0.001,coef.max,by = 0.01)
  pheatmap(coefs,cluster_rows = cluster_rows, cluster_cols = cluster_cols, show_colnames = T,col = c(colorRampPalette(c("#053061", "white"))(length(breaksList1)),colorRampPalette(c( "white","firebrick3"))(length(breaksList2))),breaks = c(breaksList1, breaksList2),display_numbers = pvls,border_color=border.col,fontsize = 10,fontsize_number = 14,number_color = "black",angle_col= angle_col,fontsize_row = 12, fontsize_col = 12)
}else{
  pheatmap(coefs,cluster_rows = cluster_rows, cluster_cols = cluster_cols, show_colnames = T,col = rev(brewer.pal(11, "RdBu")),display_numbers = pvls,border_color=border.col,fontsize = 10,fontsize_number = 14,number_color = "black",angle_col= angle_col,fontsize_row = 12, fontsize_col = 12,breaks =seq(-max,max,length.out = 12))
}
  #hm = pheatmap(coefs, cluster_cols = F, show_colnames = T, color =colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(my.breaks )),breaks = my.breaks)
  dev.copy2pdf(file = out.file, width = width, height= height,useDingbats =F)
  #coefs
}