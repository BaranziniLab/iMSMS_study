diffPlot = function(species.diff= species.difference, meta =seqmeta[seqmeta$disease !="Control", ], group = "disease",signspecies =signspecies, plot.type = c("heatmap", "barplot"), color = disease.cols){
  
  signs= unique(unlist(signspecies))
  meta[,group] = as.character(meta[,group])
  diff = species.diff[rownames(species.diff) %in% signs, match(meta$household, colnames(species.diff))]
  
  rownames(diff) = gsub("*.*s__", "", rownames(diff))
  
  if(plot.type =="heatmap"){
    annotation = data.frame(site = factor(meta$site))
    rownames(annotation) = meta$household
    anno_cols = list(site = Site.cols)
    if(!is.null(group)){
      annotation[,group] = factor(meta[,group])
      anno_cols[[group]] = color
      meta = meta[order(meta[,group]),]
    }
    
    my.breaks <- c(seq(min(species.diff), max(species.diff), by=0.01))
    
    hm = pheatmap(diff, annotation_col = annotation, annotation_colors = anno_cols, cluster_cols = F, show_colnames = F, color =colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(my.breaks )),breaks = my.breaks)
  }else if (plot.type == "barplot"){
    meta2 = split(meta, meta[,group])
    average.diff = lapply(meta2, function(x){
      data = diff[,colnames(diff) %in% x$household]
      avg = apply(data, 1, mean)
      sem = apply(data, 1, function(x) {sd(x)/sqrt(length(x))})
      smr = data.frame(mean = avg, sem =sem, group = unique(x[,group]),stringsAsFactors = F)
      rownames(smr) = rownames(data)
      smr = smr[rownames(smr) %in% signspecies[[unique(x[,group])]], ]
      smr$species = rownames(smr)
      smr
    })
    average.diffs=do.call("rbind", average.diff)
    average.diffs = average.diffs[order(average.diffs$mean), ]
    
    average.diffs$species = gsub("*.*s__", "", average.diffs$species)
    average.diffs$species = factor(average.diffs$species, level =unique(average.diffs$species))
    
    ggplot(average.diffs, aes(x=species, y= mean, fill = group )) +
      geom_bar(stat='identity',position="dodge") +
      geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                    position=position_dodge(.9))+
      coord_flip() + ylab("Difference of arcsin(relative abundance)") + xlab("")+
      scale_fill_manual(name =NULL, values = color)
    
  }
}
