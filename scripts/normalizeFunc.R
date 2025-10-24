# normalize the microbiome compostion s
# several options
# VST: Variance Stabilizing Transformation (DESeq2)
# TMM:  (edgeR)
# ILR : Isometric log ratio transform, add pseudocount 1 (compositions), n columns transformed to n-1 columns (However the interpretation of the results may be difficult, since there is no one-to-one relation between the original parts and the transformed variables.)
# CLR: center logratio transform,add pseudocount 1 (compositions)
# AST: arcsin transformation of relative abundance

normalizeFunc = function(data, method = c("rarefy", "relative", "VST","TMM", "CLR", "AST"), 
                         seqmeta = seqmeta1152){
  methods = c("rarefy", "relative", "VST","TMM", "ILR", "CLR", "AST")
  method = match.arg(method, methods)
  
  if(method == "rarefy"){
    require("GUniFrac")
    sums = colSums(data)
    qt = quantile(sums,seq(0,1, 0.01))
    depth = floor(as.numeric(qt[2]))
    data = Rarefy(data, depth = depth)
    data = data$otu.tab.rff
  }
  data.rel = sweep(data, 2, colSums(data), `/`)
  
  if(method == "relative"){
    data = data.rel
  }
  if(method == "AST"){
    data = apply(data.rel, c(1,2), function(x) {asin(sqrt(x))})
  }
  if(method  == "VST"){
    require(DESeq2)
    data = DESeqDataSetFromMatrix(countData =  data,
                                  colData = seqmeta,
                                  design = ~ disease)
    rmean =rowMeans(counts(data, normalized=F))
    len = length(rmean[rmean >=5])
    data = vst(data, nsub = len)
    data = assay(data)
  }
  
  if(method %in% c("TMM", "RLE")){
    require("edgeR")
    expr = DGEList(counts = data, genes = rownames(data))
    if(method == "TMM"){
      expr.norm = calcNormFactors(expr, method="TMM")
      data =cpm(expr.norm,normalized.lib.sizes=TRUE)
    }else{
      expr.norm = calcNormFactors(expr, method="RLE")
      data =cpm(expr.norm,normalized.lib.sizes=TRUE)
    }
  }
  
  require("compositions")
  if(method == "CLR"){
    data.tf = clr(data + 1)
    colnames(data.tf) = colnames(data)
    data = data.tf
  }
  if(method == "ILR"){
    data.tf = ilr(data + 1)
    colnames(data.tf) = colnames(data)
    data = data.tf
  }
  data
}