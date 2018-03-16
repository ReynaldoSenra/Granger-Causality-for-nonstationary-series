# Requires the packages vars and aod

toda_yamamoto<-function(x, n = 5, w = 1, lag.max = 10, type = c("const", "trend", "both", "none"), m = 1, o = 7, corrhet = 3, season = NULL, exogen = NULL, digits = 4, probability = 0.05) { 
  tablavar<-matrix(nrow = (n-1)*(ncol(x)-n-w+1), ncol = (10+2*w+2*corrhet)) 
  tablanames<-matrix(nrow = (n-1)*(ncol(x)-n-w+1), ncol = (2+w))
  cusumvar<-list()
  rootvar<-list()
  wald1i<-c()
  wald2i<-c()
  nombre<-c()
  godfrey<-c()
  white<-c()
  nomine<-colnames(x)
  datow<-matrix(nrow = nrow(x), ncol = w)
  for (i in (n+1):(n+w)){
    nombre[[i-n]]<-nomine[i]
    datow[,i-n]<-x[,i]
  }
  for (i in 2:n){
    datavar<-cbind(x[,i], datow)
    optlag<-VARselect(datavar, lag.max, type, season, exogen)
    if (nrow(x) <= 59) {
      ellag<-optlag[[1]][1]
      firstlag<-optlag[[1]][1]
      } else if (nrow(x) >= 60 & nrow(x) < 120) {
      ellag<-optlag[[1]][3]
      firstlag<-optlag[[1]][3]
      } else {
      ellag<-optlag[[1]][2]
      firstlag<-optlag[[1]][2]
    }
    datavar<-as.data.frame(datavar)
    elvar<-VAR(datavar, p = ellag, type, exogen)
    for (s in 1:corrhet){
    breush<-serial.test(elvar, lags.bg = s, type = "ES")
    hetero<-arch.test(elvar, lags.single = 0, lags.multi = s, multivariate.only = TRUE)
    godfrey[s]<-round(breush$serial[[3]], digits)
    white[s]<-round(hetero$arch.mul$p.value[[1]], digits)
    }
    if (any(godfrey <= probability) | any(white <= probability)) {
      repeat {
        ellag<-ellag+1
        elvar<-VAR(datavar, p = ellag, type, exogen)
        for (s in 1:corrhet){
          breush<-serial.test(elvar, lags.bg = s, type = "ES")
          hetero<-arch.test(elvar, lags.single = 0, lags.multi = s, multivariate.only = TRUE)
          godfrey[s]<-round(breush$serial[[3]], digits)
          white[s]<-round(hetero$arch.mul$p.value[[1]], digits)
        }
         if ((ellag == (firstlag+o)) | (any(godfrey > probability)  & any(white > probability))){ 
            if (any(godfrey <= probability) | any(white <= probability)) {
            warning(paste("Significant Autocorrelation and/ or Heteroskedasticity problems in the initial VAR of", nomine[i], "with", nombre))
          }
          break
        }
      }
    }
    elvar<-VAR(datavar, p = ellag, type, exogen)
    for (s in 1:corrhet){
      breush<-serial.test(elvar, lags.bg = s, type = "ES")
      hetero<-arch.test(elvar, lags.single = 0, lags.multi = s, multivariate.only = TRUE)
      godfrey[s]<-round(breush$serial[[3]], digits)
      white[s]<-round(hetero$arch.mul$p.value[[1]], digits)
    }
    bera<-normality.test(elvar, multivariate.only = TRUE)
    cusumvar[[(i-2)*(ncol(x)-n-w+1)+1]]<-stability(elvar, type = "OLS-CUSUM", h = 0.15, dynamic = FALSE, rescale = TRUE)
    rootvar[[(i-2)*(ncol(x)-n-w+1)+1]]<-roots(elvar, modulus = TRUE)  
    var_TD<-VAR(datavar, p = (ellag + m), type, exogen) 
    tablavar[(i-2)*(ncol(x)-n-w+1)+1,1:(8+2*corrhet)]<-c(ncol(datavar), optlag[[1]][1], optlag[[1]][2], optlag[[1]][3], optlag[[1]][4], firstlag, ellag, godfrey, white, round(bera$jb.mul$Kurtosis$p.value, digits))
    for (s in 1:w){
      wald1<-wald.test(b=coef(var_TD$varresult[[1]]), Sigma=vcov(var_TD$varresult[[1]]), Terms=seq((s+1),((w+1)*ellag),(w+1)))
      wald2<-wald.test(b=coef(var_TD$varresult[[(s+1)]]), Sigma=vcov(var_TD$varresult[[(s+1)]]), Terms=seq(1,((w+1)*ellag),(w+1)))
      tablavar[(i-2)*(ncol(x)-n-w+1)+1,8+2*corrhet+s]<-round(wald1$result[[1]][3], digits)
      tablavar[(i-2)*(ncol(x)-n-w+1)+1,8+2*corrhet+w+s+1]<-round(wald2$result[[1]][3], digits)
    }          
    tablanames[(i-2)*(ncol(x)-n-w+1)+1,1:(2+w)]<-c(nomine[i], nomine[(n+1):(n+w)], NA)
    for (j in (n+w+1):ncol(x)){
      datavar<-cbind(x[,i], datow, x[,j])
      optlag<-VARselect(datavar, lag.max, type, season, exogen)
      if (nrow(x) <= 59) {
        ellag<-optlag[[1]][1]
        firstlag<-optlag[[1]][1]
        }  else if (nrow(x) >= 60 & nrow(x) < 120) {
        ellag<-optlag[[1]][3]
        firstlag<-optlag[[1]][3]
        } else {
        ellag<-optlag[[1]][2]
        firstlag<-optlag[[1]][2]
      }
      datavar<-as.data.frame(datavar)
      elvar<-VAR(datavar, p = ellag, type, exogen)
      for (s in 1:corrhet){
        breush<-serial.test(elvar, lags.bg = s, type = "ES")
        hetero<-arch.test(elvar, lags.single = 0, lags.multi = s, multivariate.only = TRUE)
        godfrey[s]<-round(breush$serial[[3]], digits)
        white[s]<-round(hetero$arch.mul$p.value[[1]], digits)
      }
      if (any(godfrey <= probability) | any(white <= probability)) {
        repeat {
          ellag<-ellag+1
          elvar<-VAR(datavar, p = ellag, type, exogen)
          for (s in 1:corrhet){
            breush<-serial.test(elvar, lags.bg = s, type = "ES")
            hetero<-arch.test(elvar, lags.single = 0, lags.multi = s, multivariate.only = TRUE)
            godfrey[s]<-round(breush$serial[[3]], digits)
            white[s]<-round(hetero$arch.mul$p.value[[1]], digits)
          }
          if ((ellag == (firstlag+o)) | (any(godfrey > probability)  & any(white > probability))){ 
            if (any(godfrey <= probability) | any(white <= probability)) {
              warning(paste("Significant Autocorrelation and/ or Heteroskedasticity problems in the initial VAR of", nomine[i], "with", nombre, "and", nomine[j]))
            }           
            break
          }
        }
      }
      elvar<-VAR(datavar, p = ellag, type, exogen)
      for (s in 1:corrhet){
        breush<-serial.test(elvar, lags.bg = s, type = "ES")
        hetero<-arch.test(elvar, lags.single = 0, lags.multi = s, multivariate.only = TRUE)
        godfrey[s]<-round(breush$serial[[3]], digits)
        white[s]<-round(hetero$arch.mul$p.value[[1]], digits)
      }
      bera<-normality.test(elvar, multivariate.only = TRUE)
      cusumvar[[(i-2)*(ncol(x)-n-w+1)+j-n-w+1]]<-stability(elvar, type = "OLS-CUSUM", h = 0.15, dynamic = FALSE, rescale = TRUE)
      rootvar[[(i-2)*(ncol(x)-n-w+1)+j-n-w+1]]<-roots(elvar, modulus = TRUE)
      var_TD<-VAR(datavar, p = (ellag + m), type, exogen)
      tablavar[(i-2)*(ncol(x)-n-w+1)+j-n-w+1,1:(8+2*corrhet)]<-c(ncol(datavar), optlag[[1]][1], optlag[[1]][2], optlag[[1]][3], optlag[[1]][4], firstlag, ellag, godfrey, white, round(bera$jb.mul$Kurtosis$p.value, digits))
      for (s in 1:(w+1)){
        wald1<-wald.test(b=coef(var_TD$varresult[[1]]), Sigma=vcov(var_TD$varresult[[1]]), Terms=seq((s+1),((w+2)*ellag),(w+2)))
        wald2<-wald.test(b=coef(var_TD$varresult[[(s+1)]]), Sigma=vcov(var_TD$varresult[[(s+1)]]), Terms=seq(1,((w+2)*ellag),(w+2)))
        tablavar[(i-2)*(ncol(x)-n-w+1)+j-n-w+1,8+2*corrhet+s]<-round(wald1$result[[1]][3], digits)
        tablavar[(i-2)*(ncol(x)-n-w+1)+j-n-w+1,8+2*corrhet+w+s+1]<-round(wald2$result[[1]][3], digits)
      }
      tablanames[(i-2)*(ncol(x)-n-w+1)+j-n-w+1,1:(2+w)]<-c(nomine[i], nomine[(n+1):(n+w)], nomine[j])
    }
  }
  tablavar<-cbind(tablanames, tablavar)
  result_TD<-as.data.frame(tablavar, optional = FALSE)
  write.csv(result_TD, "result_Toda_Yamamoto.csv")
  names(result_TD)<-c("'Main' variable", rep("One of the variables included in w", times = w), "Additional variable", "Number of variables involved", "Akaike Information Criterion", "Hannan Quinn", "Schwarz Criterion", "FPE", "Lag indicated by the best criterion (considering the sample size)", "Lag order that corrects residual autocorrelation and Heteroskedasticity", rep("Breush Godfrey LM test for autocorrelation (p-value)", times = corrhet), rep("ARCH Heteroscedasticity test (p-value)", times = corrhet), "Jarque Bera Normality test (p-value)", rep("Does this one of the variables included in w cause the 'main' variable", times = w), "Does the 'additional' variable cause the 'main' variable", rep("Does the 'main' variable cause this one of the variables included in w", times = w), "Does the 'main' variable cause the 'additional' variable")
  View(result_TD)                                                                                                                                                                                                                                                                                                                                                                                                     
  for (r in 1:((n-1)*(ncol(x)-n-w+1))){ 
    plot(rootvar[[r]])
  }
  for (t in 1:((n-1)*(ncol(x)-n-w+1))){
    plot(cusumvar[[r]])
  }
  getwd()
}