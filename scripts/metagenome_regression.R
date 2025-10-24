#!/usr/bin/env Rscript

# To run: Rscript metagenome_regression.R humann2_metaphlan_bugs_list_256pairs_maaslin_MS_all.txt  10 disease site NULL humann2_metaphlan_bugs_list_256pairs_disease_nogenederadj.txt
args <- commandArgs(TRUE)
input=toString(args[1])
startcol=as.numeric(args[2])
type=toString(args[3])
covariate=toString(args[4])
normalizevar=toString(args[5])
output=toString(args[6])


# Format numbers function
formatnum<-function(x){
  return(formatC(signif(x,digits=3), digits=3,format="fg", flag="#"))
}

formatsci<-function(x){
  return(format(x,scientific=TRUE,digits=2))
}

roundnum<-function(x){
  return(round(x,3))
}

data=read.table(input,header=T,sep="\t")
rnames <- data[,1]
variables <- data.matrix(data[,startcol:ncol(data)])

  controls = data.matrix(data[data[,type] == "Control" ,startcol:ncol(data)])
  cases = data.matrix(data[data[, type] !="Control",startcol:ncol(data)])

logvariables=log10(variables+1e-10)
rownames(variables) <- rnames
numVARIABLEs=length(variables[1,])

# Get covariates and phenotype
names = paste("pair", 1:(nrow(data)/2),sep = "")
for(i in 1:(nrow(data)/2)){
  pair= (data$type == i)*1
  assign(names[i],pair)
}

pheno=(data[,type] !="Control")*1

# Get normalization factor (divide by AGS, multiply by geome equivalent, or multiply by 1)
if (normalizevar=="AvgGenomeSize"){
  normvar=1/data$AvgGenomeSize
  casesnorm=1/data[data[,type] !="Control",]$AvgGenomeSize
  controlsnorm=1/data[data$Pheno =="Control",]$AvgGenomeSize
}else if (normalizevar=="GenomeEquiv"){
  normvar=data$GenomeEquiv
  casesnorm=data[data[, type] !="Control",]$GenomeEquiv
  controlsnorm=data[data[, type] =="Control",]$GenomeEquiv
}else{
  normvar=rep(1,length(data[,1]))
  casesnorm=rep(1,sum(data[,type]!="Control"))
  controlsnorm=rep(1,sum(data[,type] =="Control"))
}

# Calculate stats and write to file
result=matrix(NA,nrow=numVARIABLEs,ncol=17)
count=1
for (j in 1:numVARIABLEs){
  variablename=colnames(variables)[j];
  variablevar=var(variables[,j])
  m=mean(variables[,j])
  
  # Linear regression on unnormalized variable
  if(covariate !="NULL"){
    if(covariate == "sex"){
      sex = (data$sex == "Male")*1
      s1=glm(as.formula(paste("pheno ~ variables[,j]", paste(names, "",collapse ="+"), "sex", sep="+")))
    }else if (covariate == "site"){
      site = data$Site
      sitenames = unique(site)
      for(name in sitenames){
        site = (data$Site == name)*1
        name = gsub(" ", "", name)
        assign(name,site)
      }
      sitenames = gsub(" ", "", sitenames)
      s1 = glm(as.formula(paste("pheno ~ variables[,j]", paste(names, "",collapse ="+"), paste(sitenames, "", collapse = "+"), sep="+")))
  }
    }else{
    s1=glm(as.formula(paste("pheno ~ variables[,j]", paste(names, "",collapse ="+"), sep="+"))); 
  }
  p1=coef(summary(s1))[2,4]; b1=coef(summary(s1))[2,1];
  suppressMessages(c1<-confint(s1))
  c1L=c1[2,1];c1H=c1[2,2]
  Mcases=mean(cases[,j])
  SDcases=sd(cases[,j])
  Mcontrols=mean(controls[,j])
  SDcontrols=sd(controls[,j])
  
  x=variables[,j]*normvar
  # Linear regression on normalized variable
  if(covariate !="NULL"){
    if(covariate == "sex"){
      sex = (data$sex == "Male")*1
      s2=glm(as.formula(paste("pheno ~ x", paste(names, "",collapse ="+"),"sex", sep="+")))
    }else if (covariate == "site"){
      site = data$Site
      sitenames = unique(site)
      for(name in sitenames){
        site = (data$Site == name)*1
        name = gsub(" ", "", name)
        assign(name,site)
      }
      sitenames = gsub(" ", "", sitenames)
      s1 = glm(as.formula(paste("pheno ~ variables[,j]", paste(names, "",collapse ="+"), "sex", paste(sitenames, "", collapse = "+"), sep="+")))
    }
    }else{
    s2=glm(as.formula(paste("pheno ~ x", paste(names, "",collapse ="+"), sep="+")))
  }
   p2=coef(summary(s2))[2,4]; b2=coef(summary(s2))[2,1];
  suppressMessages(c2<-confint(s2))
  c2L=c2[2,1];c2H=c2[2,2]
  McasesN=mean(cases[,j]*casesnorm)
  SDcasesN=sd(cases[,j]*casesnorm)
  McontrolsN=mean(controls[,j]*controlsnorm)
  SDcontrolsN=sd(controls[,j]*controlsnorm)
  
  result[count,]=c(formatsci(p1),formatnum(b1),formatnum(c1L),formatnum(c1H),formatsci(Mcases),formatsci(SDcases),formatsci(Mcontrols),formatsci(SDcontrols),formatsci(p2),formatnum(b2),formatnum(c2L),formatnum(c2H),formatsci(McasesN),formatsci(SDcasesN),formatsci(McontrolsN),formatsci(SDcontrolsN),variablename)
  #result[count,]=c(formatsci(p1),formatnum(b1),formatnum(c1L),formatnum(c1H),formatsci(Mcases),formatsci(SDcases),formatsci(Mcontrols),formatsci(SDcontrols),variablename)
  count=count+1;
  }

# Write results to file
colnames(result)=c("P","B","B_L","B_H","M_case","SD_case","M_cont","SD_cont","P_norm","B_norm","B_L_norm","B_H_norm","M_cases_norm","SD_cases_norm","M_cont_norm","SD_cont_norm","VARIABLE")
#colnames(result)=c("P","B","B_L","B_H","M_case","SD_case","M_cont","SD_cont","VARIABLE")
write.table(result,file=output,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")



