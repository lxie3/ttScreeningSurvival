ttscreening_surv=function(Survivaldata,DNAm,Covariate_con=NULL,Covariate_cat=NULL, Cutoff.Joint=0.1,Iterations=100,Train.SigLevel=0.05,Test.SigLevel=0.05,Percentage=0.50){

  library(fastDummies)
  library("survival")

  prepare_cat<-function(Covariate_cat){
    no.cat<-as.numeric(ncol(Covariate_cat))
    for (i in 1:no.cat){
      Covariate_cat<-dummy_cols(Covariate_cat, select_columns = colnames(Covariate_cat)[i])
      Covariate_cat<-Covariate_cat[,1:ncol(Covariate_cat)-1]
    }
    Covariate_cat<-Covariate_cat[,-(1:no.cat)]
    return(Covariate_cat)
  }

  formula_coxph<-function(alldata){
    if ( ncol(alldata)==ncol(DNAm) ) {
      formula<-as.formula(paste0("Surv(time, status)~DNAmcpg"))   #no covariate only DNAm
    } else {
      Covariate<-alldata[,-(col(DNAm))]
      formula<-as.formula(paste0("Surv(time, status)~DNAmcpg+",paste(colnames(Covariate), collapse = "+")))
    }
    return(formula)
  }

  #p value using coxph for survivaldata
  surv_coxph_p=function(datanew,formula){

    pvalue<-summary(coxph(formula, data=datanew))$coef[1,c("Pr(>|z|)")]
    return(pvalue)

  }

  #prepare data function, this data function return all the covariates in the regression including the DNAm, and DNAm is the first column
  if ( is.null(Covariate_con) & is.null(Covariate_cat)) {
    dataset<-as.data.frame(DNAm)
  } else if ( is.null(Covariate_cat)) {
    dataset<-as.data.frame(cbind(DNAm,Covariate_con))
  } else if ( is.null(Covariate_con)) {
    Covariate_cat_dummy<-prepare_cat(Covariate_cat)
    dataset<-as.data.frame(cbind(DNAm,Covariate_cat_dummy))
  } else {
    Covariate_cat_dummy<-prepare_cat(Covariate_cat)
    dataset<-as.data.frame(cbind(DNAm,Covariate_con,Covariate_cat_dummy))
  } #this step we have got all the DNAm and covariates in one dataset
  #print(ncol(dataset))
  formula<-formula_coxph(dataset)
  #print(formula)#this step we have got all the DNAm and covariates in one dataset


  p=dim(DNAm)[2]#number of potential mediators
  n=dim(DNAm)[1]#sample size
  joint.out.all=rep(NA, p)

  ### step 1 :screening each cg
  for (i in 1:p){
    tempdata<-as.data.frame(cbind(Survivaldata,DNAm[,i],dataset[,-(1:p)]))
    colnames(tempdata)[1:3]=c("time","status","DNAmcpg")
    joint.out.all[i]=surv_coxph_p(tempdata,formula)
  }

  #step 2: use 0.1 as cutoff
  Final_joint=matrix(NA,nrow=p,ncol=2)
  colnames(Final_joint)<-c("pass","tt_rate")
  rownames(Final_joint)<-colnames(DNAm)

  for (i in 1:p){
    if (joint.out.all[i] >=Cutoff.Joint) {
      Final_joint[i]=0
    }else{
      tempdata<-as.data.frame(cbind(Survivaldata,DNAm[,i],dataset[,-(1:p)]))
      colnames(tempdata)[1:3]=c("time","status","DNAmcpg")

      sig_path=rep(NA,Iterations)

      for (k in 1:Iterations) {
        set.seed(k)

        sample <- sample.int(n = n, size = floor(.67*n), replace = F)## Every loop generate IDs for train and test datasets

        ### we need to try each
        train <- tempdata[sample, ]
        test  <- tempdata[-sample, ]

        Path.out.train=surv_coxph_p(train,formula)
        Path.out.test=surv_coxph_p(test,formula)


        sig_path[k] <- ifelse (Path.out.train<Train.SigLevel &Path.out.test<Test.SigLevel,1,0)
      }
      #Final_joint[i,1]=ifelse (sum(sig_path)>Percentage*Iterations,1,0) # this Percentage*Iterations is a percentage of passing rate,default value of Percentage*Iterations = 0.5, 0<q<=1,
      Final_joint[i,1]=ifelse (sum(sig_path)>Percentage*Iterations,1,0)
      Final_joint[i,2]=sum(sig_path)/Iterations
    }
  }


  selected<-subset(Final_joint,Final_joint[,1]==1)
  No.selected<-dim(selected)[1] #number of identified mediators after ttscreening


  if (No.selected==0){
    print("No significant cpgs identified from the ttscreening method")
  }else{

    selected_ID=rownames(selected)

    DNAm_selected<-DNAm[,selected_ID]


    #two datasets: Survivaldata, DNAm_seleted

    Results<-matrix(NA,nrow=No.selected,ncol=5)
    rownames(Results)<-rownames(selected)
    colnames(Results)<-c("coef est","exp(coef)", "se(coef)", "z","Pr(>|z|)")
    for (i in 1:No.selected){
      tempdata<-as.data.frame(cbind(Survivaldata,DNAm[,i],dataset[,-(1:p)]))
      colnames(tempdata)[1:3]=c("time","status","DNAmcpg")
      #sfit <- coxph(Surv(time, status)~DNAmcpg, data=tempdata)
      if(ncol(dataset)==ncol(DNAm)){
        sfit <- coxph(as.formula(paste0("Surv(time, status)~DNAmcpg")), data=tempdata)
        Results[i,1:5]<-coef(summary(sfit))["DNAmcpg",]}
      else{sfit <- coxph(as.formula(paste0("Surv(time, status)~",paste(colnames(tempdata[,-(1:2)]), collapse = "+"))), data=tempdata)
      Results[i,1:5]<-coef(summary(sfit))["DNAmcpg",]}
    }
    Final_results<-cbind(selected,Results) #this is the results table we want the user to get eventually
    Final_results
    return(Final_results)
  }
}

