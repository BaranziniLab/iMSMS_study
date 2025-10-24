#' Plot microbial alpha-diversity
#' @diversity microbial diversity measured by shannon index and chao index
#' @meta metadata of samples
#' @group variable to shown on x-axis

diversity_plot = function(diversity, meta, group = "disease",levels = group.levels, color = group.colors,plot.type = c("barplot", "boxplot"), 
                             angle = 0, las =1, hjust =0.5, width =6, height =5,out.dir,cex = 1,jitter = 0.1){

  if(colnames(meta)[1] != "SampleID"){
    colnames(meta)[1] ="SampleID"
  }

  diversity2 = list(merge(diversity, meta, by = "SampleID"))
  
  if (plot.type == "boxplot"){
    pdf(paste(out.dir, "Alpha_diversity_", group,"_",plot.type, ".pdf",sep=""), width =width, height =height, useDingbats =F)
    layout(matrix(1:2,nrow =1))
    if(!is.null(levels)){
      if(is.list(levels)){
        diversity2[[1]][,group] = factor(diversity2[[1]][,group], levels = levels[[group]])
      }else{
        diversity2[[1]][,group] = factor(diversity2[[1]][,group], levels = levels)
      }
    }
    if(!is.null(color)){
      if(is.list(color)){
        col = color[[group]]
      }else{
        col = color
      }
      # ggplot(diversity2[[1]], aes_string(x=group, y="shannon", color="site")) +
      #   geom_boxplot() + 
      #   geom_jitter(position=position_jitter(0.2))
      
       
      boxplot(as.formula(paste("shannon~", group,sep="")), data = diversity2[[1]],
              xlab = "",ylab = "Shannon index",border= "black",las =las, boxwex = 0.4, outline =F, ylim=c(0,7))
      stripchart(as.formula(paste("shannon~", group,sep="")), data = diversity2[[1]],method = "jitter",
                 jitter = jitter, pch = 16, vertical =T,cex = cex, add=T,col = col)
      
       boxplot(as.formula(paste("chao1~", group,sep="")), data = diversity2[[1]],
             xlab = "",ylab = "Chao1 index",cex.lab = 1.3,border = "black",las = las,boxwex = 0.4, outline =F,ylim=c(0,400))
       stripchart(as.formula(paste("chao1~", group,sep="")), data = diversity2[[1]],method = "jitter",
                  jitter = jitter, pch = 16, vertical =T,cex = cex, add=T,col = col)
       
    }else{
      boxplot(as.formula(paste("shannon~", group,sep="")), data = diversity2[[1]], 
              xlab = "",ylab = "Shannon index",border= "black",las =las, boxwex = 0.4, outline =T,ylim=c(0,7))
      stripchart(as.formula(paste("shannon~", group,sep="")), data = diversity2[[1]],method = "jitter",
                 jitter = jitter, pch = 16, vertical =T,cex =cex, add=T)
      
      boxplot(as.formula(paste("chao1~", group,sep="")), data = diversity2[[1]],
              xlab = "",ylab = "Chao1 index",cex.lab = 1.3,border = "black",las = las,boxwex = 0.4, outline =T,ylim=c(0,400))
      stripchart(as.formula(paste("chao1~", group,sep="")), data = diversity2[[1]],method = "jitter",
                 jitter = jitter, pch = 16, vertical =T,cex =cex, add=T )
    }
      
    dev.off()
  }else if(plot.type == "barplot"){
    data_summary = function(data, varname, groupnames){
      require(plyr)
      summary_func = function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE),
          sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]][!is.na(x[[col]])])))
      }
      data_sum = ddply(data, groupnames, .fun=summary_func,
                       varname)
      data_sum = rename(data_sum, c("mean" = varname))
      return(data_sum)
    }
    shannon_sum = lapply(diversity2, function(x){
      y = x[,c(group,"shannon")]
      y = data_summary(y, varname="shannon", groupnames =group) 
      y
    })
    
    chao_sum = lapply(diversity2,function(x){
      y = x[,c(group,"chao1")]
      y =data_summary(y, varname="chao1", groupnames = group)
    } )
    
    colnames(shannon_sum[[1]])[1] = "disease" 
    colnames(chao_sum[[1]])[1] = "disease" 
    
    ### plot shannon and chao1 in one figure
    shannon_sum2 = lapply(shannon_sum, function(x){
      x$Description = "Shannon"
      x
    })
    chao_sum2 = lapply(chao_sum, function(x){
      x$Description = "Chao1"
      x
    })
    
    if(!is.null(levels)){
      shannon_sum2[[1]]$disease = factor(shannon_sum2[[1]]$disease, levels = levels[[group]])
      chao_sum2[[1]]$disease = factor(chao_sum2[[1]]$disease, levels = levels[[group]])
    }
    
    if(!is.null(color)){
      if(is.list(color)){
        col = color[[group]]
      }else{
        col = color
      }
      p1 = ggplot(shannon_sum2[[1]],aes(x=disease,y= shannon, fill =disease)) + 
        geom_bar(stat="identity", position=position_dodge(), width= 0.5,show.legend = F,color="black") +
        geom_errorbar(aes(ymin=shannon, ymax =shannon + sem), width=.2,
                      position=position_dodge(.5))+ 
        scale_fill_manual(values = col)+
        xlab("")+ ylab("Shannon index")+
        theme(axis.text.x = element_text(angle= angle,hjust =hjust))
      
      p2 = ggplot(chao_sum2[[1]],aes(x=disease, y = chao1, fill =disease)) + 
        geom_bar(stat="identity", position=position_dodge(), width= 0.5,show.legend = F,color="black") +
        geom_errorbar(aes(ymin=chao1, ymax =chao1 + sem), width=.2,
                      position=position_dodge(.5)) + 
        scale_fill_manual(values = col)+
        xlab("")+ylab("Chao1 index")+
        theme(axis.text.x = element_text(angle=angle,hjust =hjust))
    }else{
      p1 = ggplot(shannon_sum2[[1]],aes(x=disease,y= shannon, fill =disease)) + 
        geom_bar(stat="identity", position=position_dodge(), width= 0.5,show.legend = F,color="black") +
        geom_errorbar(aes(ymin=shannon, ymax =shannon + sem), width=.2,
                      position=position_dodge(.5))+ 
        xlab("")+ ylab("Shannon index")+
        theme(axis.text.x = element_text(angle= angle,hjust =hjust))
      
      p2 = ggplot(chao_sum2[[1]],aes(x=disease, y = chao1, fill =disease)) + 
        geom_bar(stat="identity", position=position_dodge(), width= 0.5,show.legend = F,color="black") +
        geom_errorbar(aes(ymin=chao1, ymax =chao1 + sem), width=.2,
                      position=position_dodge(.5)) + xlab("")+ylab("Chao1 index")+
        theme(axis.text.x = element_text(angle=angle,hjust =hjust))
    }
    
    if(!require("gridExtra")){
      install.packages("gridExtra")
    }
    require("gridExtra")
    g = grid.arrange(p1, p2, nrow = 1)
    ggsave(g, file = paste(out.dir, "Alpha_diversity_", group,"_",plot.type, ".pdf",sep=""), width =width, height =height, useDingbats =F)
    diversity2[[1]]
  } 
}