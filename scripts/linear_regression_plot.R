# @linear.res output of the linear regression result by {linear_regression} function
# @disease.only display only disease/disease course but not age or sex or bmi
# @disease.each linear regression on each diease course
# @order.group group used to order the coefficient
# @groups number of groups to plot 

require(ggplot2)
require(reshape2)
require(reshape)
linear_regression_plot =function(linear.res,disease.only =TRUE,group = "disease",taxaname=T, order.group = "MS",disease.each =FALSE, fdr = 0.05, taxa = "ASV",col.disease = "red", 
                                 out.file, padjust =TRUE, width =10,height =12, title= ""){
  if(taxaname){
    taxaname = sapply(linear.res$taxonomy,function(x){
      y = gsub("; [a-z]__;.*", "", x)
      if(grepl("s__",y)){
        y = gsub("^.*; g__", "g__",y) 
      }else{
        strs = unlist(strsplit(y, "; "))
        y = strs[length(strs)]
        if(grepl("__$", y)){
          y = strs[length(strs)-1]
        }
      }
      y
    })
    
    if(taxa  == "ASV"){
      taxaname = paste(linear.res$ASVID, ": ", taxaname, sep="")
    }else{
      taxaname = linear.res$ID
    }
    linear.res$ID = taxaname
  }else{
    linear.res$ID =linear.res$taxonomy
    
  }
  
  res = linear.res[order(linear.res[,paste("Coef_", group, order.group,sep="")],decreasing = T),]
  
  if(disease.only){
    if(group =="disease"){
      diseases = order.group
      col.selects = c(colnames(res)[grepl(diseases,colnames(res))], "ID", "taxonomy")
      col.selects = col.selects[!grepl("Pr_|fdr_", col.selects)]
      if(padjust){
        res.sign1 = res[res[,paste("fdr_disease", diseases,sep="")] < fdr ,match(col.selects, colnames(res))]
      }else{
        res.sign1 = res[res[,paste("Pr_disease", diseases,sep="")] < fdr ,match(col.selects, colnames(res))]
      }
    }else{
      if(disease.each){
        diseases = order.group
        col.selects = c(colnames(res)[grepl(paste(diseases,collapse = "|"),colnames(res))], "ID", "taxonomy")
        col.selects = col.selects[!grepl("Pr_|fdr_", col.selects)]
        if(padjust){
          res.sign1 = res[res[,paste("fdr_disease_course", order.group,sep="")] < fdr ,match(col.selects, colnames(res))]
        }else{
          res.sign1 = res[res[,paste("Pr_disease_course", order.group,sep="")] < fdr ,match(col.selects, colnames(res))]
        }
      }else{
        diseases = c("RRMS","SPMS","PPMS")
        col.selects = c(colnames(res)[grepl(paste(diseases,collapse = "|"),colnames(res))], "ID", "taxonomy")
        col.selects = col.selects[!grepl("Pr_|fdr_", col.selects)]
        if(padjust){
          res.sign1 = res[res$fdr_disease_courseRRMS < fdr | res$fdr_disease_courseSPMS < fdr | res$fdr_disease_coursePPMS < fdr ,match(col.selects, colnames(res))]
        }else{
          res.sign1 = res[res$Pr_disease_courseRRMS < fdr | res$Pr_disease_courseSPMS < fdr  | res$Pr_disease_coursePPMS < fdr ,match(col.selects, colnames(res))]
        }
      }
    }
    
    
    if(nrow(res.sign1) > 0){
      colnames(res.sign1) = c(diseases, as.vector(t(outer(diseases, c("_coef_low", "_coef_high"), paste, sep=""))), "ID", "taxonomy")
      coef = res.sign1[,c("ID", diseases)]
      coef2 = melt(coef, id = "ID")
      
      coef_low = res.sign1[,c("ID", colnames(res.sign1)[grepl("low",colnames(res.sign1))])]
      coef_low2 = melt(coef_low,id= "ID")
      colnames(coef_low2)[3] = "coef_low"
      coef_low2$variable = gsub("_coef_low", "",coef_low2$variable)
      
      coef_high = res.sign1[,c("ID", colnames(res.sign1)[grepl("high",colnames(res.sign1))])]
      coef_high2= melt(coef_high, id = "ID")
      coef_high2$variable = gsub("_coef_high", "",coef_high2$variable)
      colnames(coef_high2)[3] = "coef_high"
      
      plot.data =merge(coef2, coef_low2, by =c("ID", "variable"))
      plot.data =merge(plot.data, coef_high2, by =c("ID", "variable"))
      plot.data$ID = factor(plot.data$ID, levels= rev(res.sign1$ID))
      plot.data$variable =factor(plot.data$variable,levels =rev(diseases))
      
      # define the face for labeled significant microbes
      faces = rep("italic",nrow(res.sign1))
      #faces[res.sign2.index] = "bold.italic"
      if(taxa != "ASV"){
        labs = gsub(paste(taxa, ".*: ",sep=""), "",rev(res.sign1$ID))
      }else{
        labs = rev(res.sign1$ID)
      }
      
      max1 = max(plot.data$coef_high)
      min1 = min(plot.data$coef_low)
      if(abs(min1) > max1){
        ylim.max = abs(min1)
        ylim.min = min1
      }else{
        ylim.max = max1
        ylim.min = -max1
      }
      g = ggplot(plot.data, aes(x = ID, y =value,color = variable)) +
        geom_point(position=position_dodge(.9), aes(y = value, colour = variable)) +
        scale_color_manual(values = rev(col.disease)) + 
        geom_errorbar(aes(ymin=coef_low, ymax=coef_high), width=0.5,position = position_dodge(0.9),size =1) +
        scale_x_discrete(breaks= rev(res.sign1$ID), labels=labs) + ylim(ylim.min,ylim.max) +
        ylab("Regression coefficient") + xlab("")+ theme(legend.position = "none") +
        ggtitle(title)+
        #theme(panel.border = element_rect(size =1)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        coord_flip()
      
      ggsave(g, file = out.file,width =width, height =height, useDingbats =F)
      g
    }else{
      return("No significant bacteria!")
    }
  }else{
    if(group =="disease"){
      diseases =order.group
      col.selects = c(colnames(res)[grepl(paste(order.group, "|Female|age|bmi",sep=""),colnames(res))], "ID", "taxonomy")
      col.selects = col.selects[!grepl("Pr_|fdr_", col.selects)]
      if(padjust){
        res.sign1 = res[res[,paste("fdr_disease",diseases, sep="")] < 0.05 ,match(col.selects, colnames(res))]
      }else{
        res.sign1 = res[res[,paste("Pr_disease",diseases, sep="")] < 0.05 ,match(col.selects, colnames(res))]
      }
    }else {
      if(disease.each){
        diseases = order.group
        col.selects = c(colnames(res)[grepl(paste(c(diseases, "Female", "age", "bmi"),collapse = "|"),colnames(res))], "ID", "taxonomy")
        col.selects = col.selects[!grepl("Pr_|fdr_", col.selects)]
        if(padjust){
          res.sign1 = res[res[,paste("fdr_disease_course", order.group,sep="")] < 0.05 ,match(col.selects, colnames(res))]
        }else{
          res.sign1 = res[res[,paste("Pr_disease_course", order.group,sep="")] < 0.05 ,match(col.selects, colnames(res))]
        }
        
      }else{
        diseases = c("RRMS","SPMS","PPMS")
        col.selects = c(colnames(res)[grepl(paste(c(diseases,"Female", "age","bmi"),collapse = "|"),colnames(res))], "ID", "taxonomy")
        col.selects = col.selects[!grepl("Pr_|fdr_", col.selects)]
        if(padjust){
          res.sign1 = res[res$fdr_disease_courseRRMS < 0.05 | res$fdr_disease_courseSPMS < 0.05  | res$fdr_disease_coursePPMS < 0.05 ,match(col.selects, colnames(res))]
        }else{
          res.sign1 = res[res$Pr_disease_courseRRMS < 0.05 | res$Pr_disease_courseSPMS < 0.05  | res$Pr_disease_coursePPMS < 0.05 ,match(col.selects, colnames(res))]
        }
       
      }
    }
    
    if(nrow(res.sign1) > 0){
      colnames(res.sign1) = c(c(diseases,"Female","age","bmi"), as.vector(t(outer(c(diseases,"Female","age","bmi"), c( "_coef_low", "_coef_high"), paste, sep=""))), "ID", "taxonomy")
      coef = res.sign1[,c("ID", c(diseases, "Female", "age", "bmi"))]
      coef2 = melt(coef, id = "ID")
      
      coef_low = res.sign1[,c("ID", colnames(res.sign1)[grepl("low",colnames(res.sign1))])]
      coef_low2 = melt(coef_low,id= "ID")
      colnames(coef_low2)[3] = "coef_low"
      coef_low2$variable = gsub("_coef_low", "",coef_low2$variable)
      
      coef_high = res.sign1[,c("ID", colnames(res.sign1)[grepl("high",colnames(res.sign1))])]
      coef_high2= melt(coef_high, id = "ID")
      coef_high2$variable = gsub("_coef_high", "",coef_high2$variable)
      colnames(coef_high2)[3] = "coef_high"
      
      plot.data =merge(coef2, coef_low2, by =c("ID", "variable"))
      plot.data =merge(plot.data, coef_high2, by =c("ID", "variable"))
      plot.data$ID = factor(plot.data$ID, levels= rev(res.sign1$ID))
      plot.data$variable =factor(plot.data$variable,levels =rev(c(diseases, "Female", "age","bmi")))
      
      # define the face for labeled significant microbes
      faces = rep("italic",nrow(res.sign1))
      #faces[res.sign2.index] = "bold.italic"
      if(taxa != "ASV"){
        labs = gsub(paste(taxa, ".*: ",sep=""), "",rev(res.sign1$ID))
      }else{
        labs = rev(res.sign1$ID)
      }
    
      g = ggplot(plot.data, aes(x = ID, y =value,color = variable)) +
        geom_point(position=position_dodge(.9), aes(y = value, colour = variable)) +
        scale_color_manual(values = rev(c(col.disease,"Female" = "#018571","age" = "#5e3c99","bmi" = "orange" ))) + 
        geom_errorbar(aes(ymin=coef_low, ymax=coef_high), width=0.5,position = position_dodge(0.9),size =1) +
        scale_x_discrete(breaks= rev(res.sign1$ID), labels=labs) +
        ylab("Regression coefficient") + xlab("")+ theme(legend.position = "none") +
        #theme(panel.border = element_rect(size =1)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        ggtitle(title)+
        coord_flip()
      ggsave(g, file = out.file,width =width, height =height, useDingbats =F,limitsize = FALSE)
      g
    }else{
      return("No significant bacteria!")
    }
  }
  plot.data
}