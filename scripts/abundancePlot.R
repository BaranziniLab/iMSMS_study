abundancePlot = function(abundance, metadata=, group1 ="disease", group2 = "allergies", feature.plot ="Akkermansia.s__Akkermansia_muciniphila", arcsin =TRUE, color = NULL, ylab = "Akkermansia muciniphila", ymax = NULL){
  if(!is.null(group2)){
    metadata = metadata[metadata[,group2] != "" ,]
  }
  
  if(arcsin){
    
    abundance = asin(sqrt(abundance))
    #abundance = log10(abundance)
  }
  
  rownames(metadata) = gsub("-",".", metadata$iMSMS_ID)
  data = merge(t(abundance), metadata, by = "row.names")
  colnames(data) = gsub("[|]| ", ".",colnames(data))
  colnames(data) = gsub("[[]|[]]|[:]", "",colnames(data))
  if(is.null(ymax)){
    ymax =max(abundance)
  }else{
    ymax = ymax
  }
  if(!is.null(color)){
    if(!is.null(group2)){
      ggplot(data, aes_string(x = group1, y = feature.plot,fill = group2))+
        geom_jitter(aes_string(color =group2), position=position_jitterdodge(0.3), size = 1)+
        ylab(ylab) + 
        ylim(0, ymax)+
        geom_boxplot( outlier.color = NA, notch=T,alpha=0.2) +
        xlab("") + scale_color_manual(values = color) + scale_fill_manual(values = color)
    }else{
      ggplot(data, aes_string(x = group1, y = feature.plot,color = group1))+
        geom_jitter(position=position_jitterdodge(jitter.width = 0.5),size = 1)+
        ggtitle(ylab) +
        ylab("Arcsin(relative abundance)") + 
        ylim(0, ymax)+
        geom_boxplot(outlier.color = NA, notch=T,alpha=0.2) +
        xlab("") +scale_color_manual(values = color) +  theme(legend.position = "none") 
    }
  }else{
    if(!is.null(group2) ){
      ggplot(data, aes_string(x = group1, y = feature.plot, fill =group2)) +
        geom_jitter(aes_string(color =group2) ,position=position_jitterdodge(0.3), size = 1 )+
        geom_boxplot(outlier.color = NA, notch=T,alpha=0.2) + ylab(ylab) + xlab("") +
        ylim(0, ymax) + 
        theme(legend.position = "none") 
    }else{
      ggplot(data, aes_string(x = group1, y = feature.plot,color =group1)) +
        geom_jitter(position=position_jitterdodge(0.3), size = 1 )+
        geom_boxplot(outlier.color = NA, notch=T,alpha=0.2) + ylab(ylab) + xlab("") + theme(legend.position = "none") 
      + ylim(0, ymax)
    }
  }
}
