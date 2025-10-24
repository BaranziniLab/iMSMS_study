lapply(list("ape", "ggplot2","phyloseq", "car", "rgl"), require, character.only =TRUE)
# line, add line to HC-MS paired samples

unifracPCA  = function(file="results/3.diversity/Filtered_otu_averaged/diversity10k/bdiv_even10000/weighted_unifrac_dm.txt",  
                      condition = meta, group = "wet_dry",  color.list=color.list, adonis =TRUE,strataID = "pair",
                      levels = c("Control","RRMS","PMS"), method =c("PCoA", "3DPCoA"),shapes =c(16,16,16),
                      width =6, height=5, line =FALSE,label=FALSE, firmicutes=FALSE, color.max = 1, low = "blue", high = "red"){
  if(is.character(file)){
    weight = read.table(file, head=T,as.is=T)
  }else{
    weight =file
  }
 
  colnames(weight) = gsub("^X", "", colnames(weight))
  pcs = cmdscale(weight, k =3, eig = T, x.ret =T, list. =T)
  pcs$eig[pcs$eig < 0 ] = 0
  # calculate the % of variance explained
  percent = as.numeric(pcs$eig/sum(pcs$eig))[1:3]
  condition = condition[match(rownames(weight), condition[,1]),]
  
  if(adonis){
    #condition2 = condition[!condition$pair %in% as.numeric(names(which(table(condition$pair)==1))), ]
    condition2 = condition
    weight2 = as.matrix(weight[rownames(weight) %in% condition2[,1],colnames(weight) %in% condition2[,1]])

    condition2 =condition2[match(rownames(weight2),condition2[,1]),]
    if(is.null(strataID)){
      adonisres = vegan::adonis(as.formula(paste("weight2 ~", group,sep="")) , data =condition2)
    }else{
      adonisres = vegan::adonis(as.formula(paste("weight2 ~", group,sep="")) , data =condition2, 
                                strata = condition2[,strataID])
    }
  }
  if(is.list(color.list)){
    cols= color.list[[group]]
  }else{
    cols = color.list
  }
  if(method == "3DPCoA"){
    p= scatter3d(pcs$points[,1],-pcs$points[,2],pcs$points[,3],
                 groups=factor(condition[,group],levels =levels),ellipsoid = F, 
                 ellipsoid.alpha = 0.8,grid = FALSE, surface = FALSE,
                 axis.col = c("black", "black", "black"), 
                 xlab = paste("PC1", "(", round(percent[1]*100, digits = 2), "%", ")",sep=""), 
                 ylab = paste("PC2", " (",round(percent[2]*100, digits = 2), "%", ")",sep=""), 
                 zlab = paste("PC3", " (",round(percent[3]*100, digits = 2), "%", ")",sep=""),
                 axis.scales = FALSE, surface.col  = cols)
    rgl.postscript(paste("UniFrac_", group, "_3D_PCoA.pdf",sep=""),fmt="pdf")
  }else{
    data = data.frame(pcs$points)
    rownames(condition) = condition[,1]
    data = merge(data, condition, by ="row.names")
    if(!is.null(levels)){
      data[,group] = factor(data[,group], levels=levels)
    }else{
      data[,group] = factor(data[,group])
    }
    
    if(line){
      p = ggplot(data = data, aes_string(x= "X1", y ="X2", color =group, shape = group)) + 
        geom_line(aes(group=pair), color ="grey", linetype = 2)+
        geom_point(aes_string(shape =group)) +scale_color_manual(values = cols) + 
        scale_shape_manual(values = shapes) + 
        xlab(paste("PC1", "(", round(percent[1]*100, digits = 2), "%", ")",sep="")) +
        ylab(paste("PC2", " (",round(percent[2]*100, digits = 2), "%", ")",sep="")) 
      
    }else{
      if(firmicutes){
      colors = color.list[[group]]
      colors = colors[match(data$Row.names,names(colors))]
      data$color = colors
        p = ggplot(data = data, aes_string(x= "X1", y ="X2")) + 
          geom_point(aes(colour = color)) +scale_color_gradient(limits=c(0,color.max), low = low, high = high) + 
          xlab(paste("PC1", "(", round(percent[1]*100, digits = 2), "%", ")",sep="")) +
          ylab(paste("PC2", " (",round(percent[2]*100, digits = 2), "%", ")",sep=""))
      }else{
        p = ggplot(data = data, aes_string(x= "X1", y ="X2", color =group, shape = group)) + 
          geom_point(aes_string(shape =group)) +scale_color_manual(values = cols) + 
          scale_shape_manual(values = shapes) + 
          xlab(paste("PC1", "(", round(percent[1]*100, digits = 2), "%", ")",sep="")) +
          ylab(paste("PC2", " (",round(percent[2]*100, digits = 2), "%", ")",sep=""))
      }
    }
    if(label){
      p = p + geom_text(aes_string(label= "map$SampleID", colour = legend),size =1.5, hjust = 0, check_overlap = T, nudge_x = 0)
    }
    ggsave(p, file= paste("UniFrac_", group, "_PCoA.pdf",sep=""), 
           width =width, height =height, useDingbats=F)
  }
  if(adonis){
    return(list(p,adonisres))
  }
  return(p)
}