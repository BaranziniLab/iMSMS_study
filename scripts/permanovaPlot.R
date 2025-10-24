permanovaPlot =function(data, confounder.anova, distance =TRUE, termgroup, category = "kegg"){
  names = colnames(confounder.anova)
  names = names[!names %in% c("iMSMS_ID", "household")]
  
  adonis.res = list()
  for(i in 1:length(names)){
    confounder2 = confounder.anova[confounder.anova$iMSMS_ID %in% rownames(data),]
    confounder2 = confounder2[!is.na(confounder2[,i+1]) & confounder2[,i+1] != "", ]
    
    if(! distance){
      abun = data[rownames(data) %in% confounder2$iMSMS_ID, ]
      dis = vegdist(abun, method = "bray",diag =T, upper = T)
    }else{
      dis = data[match(confounder2$iMSMS_ID, rownames(data)), match(confounder2$iMSMS_ID, colnames(data))]
    }
    adonis.res[[i]] = vegan::adonis(as.formula(paste("dis ~",names[i], sep = "")), data = confounder2)
  }
  names(adonis.res) = names
  # extract the R2 and Pvalue
  extra.num = 0 # number of above house site covariate
  
  result = matrix(NA, nrow = length(names)+extra.num, ncol =2)
  for(i in 1:(length(names)+extra.num)){
    result[i,1] = adonis.res[[i]]$aov.tab$R2[1]
    result[i,2] = adonis.res[[i]]$aov.tab$`Pr(>F)`[1]
  }
  
  rownames(result) = c(names)
  colnames(result) = c("R2", "Pvalue")
  result = data.frame(result, stringsAsFactors = F)
  result$Padjust = p.adjust(result$Pvalue, method = "fdr")
  
  # write.csv(result, "Adonis_R2.csv", )
  # result=read.csv("Adonis_R2.csv",head=T,as.is=T,row.names = 1)
  result$ID = rownames(result)
  for(i in 1:nrow(result)){
    result$Group[i] = termgroup[termgroup$Term == result$ID[i], "Group"]
  }
  presult = result[result$Pvalue < 0.05,]
  padj.result = result[result$Padjust < 0.05,]
  anova.cols = c("Demography" = "#0571B0", "Disease"="#CCCCCC", "Life style" ="#92C5DE", "Medication"="#F4A582", "Physiology"="#CA0020")
  g = ggplot(result, aes(x = reorder(ID, R2),y=R2, fill = Group)) +
    geom_bar(stat='identity') +
    coord_flip() + ylab("Adonis R2") + xlab("") +
    scale_fill_manual(values = anova.cols) +
    geom_text(data = presult, aes(ID, R2),label="*", col= "black",nudge_y = 0.005, nudge_x = -0.15)+
    geom_text(data = padj.result, aes(ID, R2),label="*", col= "red",nudge_y = 0.01,nudge_x = -0.15)
  
  ggsave(g, file = paste(category, "_PERMANOVA_weighted_unifrac_Adonis_R2_barplot.pdf",sep=""), width =7.5, height =5.5, useDingbats =F)
  g
}

