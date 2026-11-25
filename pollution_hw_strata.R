library(data.table)
library(lubridate)
library(splines)
library(rmeta)
library(mixmeta)
library(dplyr)
library(mgcv)

load("/Users/sunweiran/Desktop/interaction_SWR/raw_data/272city_0216.Rdata")

yg.hw<-function(v,days,q)# v is your temperature; days is heat waves days;q is quantile;
{
  cut_point<-quantile(v,q, na.rm=T)
  s.hw<-rep(0,length(v))
  s.hw[v>=cut_point]<-1
  s.hw.sum<-Lag(s.hw,c(0:(days-1))) %>% rowSums()
  s.hw.sum[s.hw.sum<days]<-0
  l.hw <- Lag(s.hw.sum,c(0:(1-days))) %>% rowSums()
  l.hw[l.hw>0]<-1
  hw<-l.hw
  hw[is.na(hw)]<-0
  return(hw)
}
#define heat wave
list_hw <- lapply(list1, function(x) x[,hw_mean_lag2_p90:= yg.hw(temp, 2, 0.9)])
list_hw <- lapply(list_hw, function(x) x[,hw_mean_lag2_p92.5:= yg.hw(temp, 2, 0.925)])
list_hw <- lapply(list_hw, function(x) x[,hw_mean_lag2_p95:= yg.hw(temp, 2, 0.95)])
list_hw <- lapply(list_hw, function(x) x[,hw_mean_lag2_p97.5:= yg.hw(temp, 2, 0.975)])
list_hw <- lapply(list_hw, function(x) x[,hw_mean_lag3_p90:= yg.hw(temp, 3, 0.9)])
list_hw <- lapply(list_hw, function(x) x[,hw_mean_lag3_p92.5:= yg.hw(temp, 3, 0.925)])
list_hw <- lapply(list_hw, function(x) x[,hw_mean_lag3_p95:= yg.hw(temp, 3, 0.95)])
list_hw <- lapply(list_hw, function(x) x[,hw_mean_lag3_p97.5:= yg.hw(temp, 3, 0.975)])
#list.cat1 <- lapply(list_hw, function(x) x[, cat.o3 := as.numeric(cut(l02o3, c(0,100,999)))])


#####O3 effcts of different levels of heatwave#####
list2 <- list()
for (d in c("lag2","lag3")){
  for (p in c("p90", "p92.5", "p95", "p97.5")) {
    hw <- paste0("hw_mean_",d,"_",p)
    list3 <- list()
    #hw <- "hw_mean_lag2_p90"
    for (j in c("total", "cvd", "res")) {
      #fml <- get(j)~l02o3+ns(l02tmp_mean,6)+ns(l02rh_mean,3)+ as.factor(dow) +ns(date,df=round(7*length(date)/365.25))
      fml <- get(j)~l02o3+ns(l02rh,3)+ as.factor(dow) +ns(date,df=round(7*length(date)/365.25))
      model1 <- list();model2 <- list()
      for(i in 1:272) {
        print(paste0(Sys.time(), " ", j, " ", i, " ", hw))
        dt <- list_hw[[i]]
        #cut <- quantile(dt$l02tmp_mean, 0.5, na.rm = T)
        tryCatch({
          model1[[i]] <- gam(fml, data = dt[get(hw) == 1,], family=quasipoisson,na.action=na.exclude)
          model2[[i]] <- gam(fml, data = dt[get(hw) == 0,], family=quasipoisson,na.action=na.exclude)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
      
      coef1<-vector()
      se1<-vector()
      var1<-vector()
      t1<-vector()
      p1<-vector()
      
      coef2<-vector()
      se2<-vector()
      var2<-vector()
      t2<-vector()
      p2<-vector()
      
      rsq<-vector()
      gcv<-vector()
      
      code_city<-vector()
      
      for (m in 1:length(model1)){
        if (!is.null(model1[[m]]) & !is.infinite(model1[[m]]$gcv.ubre)){
          coef1[m]<-summary(model1[[m]])$p.coeff[2]
          se1[m]<-summary(model1[[m]])$se[2]
          var1[m]<- se1[m]^2
          t1[m]<-summary(model1[[m]])$p.t[2]
          p1[m]<-summary(model1[[m]])$p.pv[2]
          rsq[m]<-summary(model1[[m]])$r.sq
          gcv[m]<-model1[[m]]$gcv
        }else{
          coef1[m]<-NA
          se1[m]<-NA
          var1[m]<- NA
          t1[m]<-NA
          p1[m]<-NA
          rsq[m]<-NA
          gcv[m]<-NA
        }
        if (!is.null(model2[[m]]) & !is.infinite(model2[[m]]$gcv.ubre)){
          coef2[m]<-summary(model2[[m]])$p.coeff[2]
          se2[m]<-summary(model2[[m]])$se[2]
          var2[m]<- se2[m]^2
          t2[m]<-summary(model2[[m]])$p.t[2]
          p2[m]<-summary(model2[[m]])$p.pv[2]
          rsq[m]<-summary(model2[[m]])$r.sq
          gcv[m]<-model2[[m]]$gcv
        }else{
          coef2[m]<-NA
          se2[m]<-NA
          var2[m]<- NA
          t2[m]<-NA
          p2[m]<-NA
          rsq[m]<-NA
          gcv[m]<-NA
        }
        code_city[m] <- list1[[m]]$code_city[1]
      }
      #print(paste0(Sys.time(), " model1 ", m))
      list3[[j]] <- data.table(outcome = j,strata = hw,code_city,coef1,se1,var1,t1,p1,coef2,se2,var2,t2,p2,rsq,gcv)
    }
    list2[[hw]] <- rbindlist(list3)
    rm(list3);gc()
  }
}


save(list2, file = "/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/list_o3_hw_strata.Rdata")
data.O3 <- rbindlist(list2)
fwrite(data.O3, "/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/list_o3_hw_strata.csv")

#meta
data.O3 <- fread("/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/list_o3_hw_strata.csv")
data.O3[abs(coef1) > 0.1 | se1 > 0.1,]$var1 <- NA
data.O3[abs(coef2) > 0.1 | se2 > 0.1,]$var2 <- NA

list3 <- list()
for (d in c("lag2","lag3")) {
  for (p in c("p90", "p92.5", "p95", "p97.5")) {
    list4 <- list()
    hw <- paste0("hw_mean_",d,"_",p)
    #hw <- "hw_mean_lag2_p90"
    for (i in c("total", "cvd", "res")){
      #i <- "total"
      a <- data.O3[outcome == i, ]
      a <- a[strata == hw, ]
      print(paste0(Sys.time(), " meta ", i, " ", hw))
      indnomiss <- which(!is.na(a$var1))
      y1 <- a[indnomiss,coef1]
      S1 <- a[indnomiss,var1]
      meta1 <- mixmeta(y1, S1, random=~1|code_city, data=data.city[indnomiss,], 
                       control=list(showiter=F))
      beta1<-meta1$coefficients
      se1<-sqrt(meta1$vcov)
      betaoverall1<-round((exp(beta1*10)-1)*100,2)
      betaoverallinf1<-round((exp((beta1*10-1.96*se1*10))-1)*100,2)
      betaoverallsup1<-round((exp((beta1*10+1.96*se1*10))-1)*100,2)
      
      indnomiss <- which(!is.na(a$var2))
      y2 <- a[indnomiss,coef2]
      S2 <- a[indnomiss,var2]
      meta2 <- mixmeta(y2, S2, random=~1|code_city, data=data.city[indnomiss,], 
                       control=list(showiter=F))
      beta2<-meta2$coefficients
      se2<-sqrt(meta2$vcov)
      betaoverall2<-round((exp(beta2*10)-1)*100,2)
      betaoverallinf2<-round((exp((beta2*10-1.96*se2*10))-1)*100,2)
      betaoverallsup2<-round((exp((beta2*10+1.96*se2*10))-1)*100,2)
      
      z_strata <- (beta2 - beta1) / sqrt(se1^2 + se2^2)
      p_strata <- 2 * pnorm(-abs(z_strata))
      
      list4[[i]] <- data.table(outcome = i, strata = hw,
                               n_hw = paste0(betaoverall2, " (",betaoverallinf2,",",betaoverallsup2, ")"),
                               hw    = paste0(betaoverall1, " (",betaoverallinf1,",",betaoverallsup1, ")"))
      
      list4[[i]]$p <- p_strata 
    }
    list3[[hw]] <- rbindlist(list4)
    rm(list4);gc()
  }
}
pool.O3 <- rbindlist(list3)
fwrite(pool.O3, "/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/pool_o3_hw_strata.csv")
rm(list3);gc()


#####PM25 effcts of different levels of heatwave#####
list2 <- list()
for (d in c("lag2","lag3")){
  for (p in c("p90", "p92.5", "p95", "p97.5")) {
    hw <- paste0("hw_mean_",d,"_",p)
    list3 <- list()
    #hw <- "hw_mean_lag2_p90"
    for (j in c("total", "cvd", "res")) {
      #fml <- get(j)~l02o3+ns(l02tmp_mean,6)+ns(l02rh_mean,3)+ as.factor(dow) +ns(date,df=round(7*length(date)/365.25))
      fml <- get(j)~l02pm25+ns(l02rh,3)+ as.factor(dow) +ns(date,df=round(7*length(date)/365.25))
      model1 <- list();model2 <- list()
      for(i in 1:272) {
        print(paste0(Sys.time(), " ", j, " ", i, " ", hw))
        dt <- list_hw[[i]]
        #cut <- quantile(dt$l02tmp_mean, 0.5, na.rm = T)
        tryCatch({
          model1[[i]] <- gam(fml, data = dt[get(hw) == 1,], family=quasipoisson,na.action=na.exclude)
          model2[[i]] <- gam(fml, data = dt[get(hw) == 0,], family=quasipoisson,na.action=na.exclude)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
      
      coef1<-vector()
      se1<-vector()
      var1<-vector()
      t1<-vector()
      p1<-vector()
      
      coef2<-vector()
      se2<-vector()
      var2<-vector()
      t2<-vector()
      p2<-vector()
      
      rsq<-vector()
      gcv<-vector()
      
      code_city<-vector()
      
      for (m in 1:length(model1)){
        if (!is.null(model1[[m]]) & !is.infinite(model1[[m]]$gcv.ubre)){
          coef1[m]<-summary(model1[[m]])$p.coeff[2]
          se1[m]<-summary(model1[[m]])$se[2]
          var1[m]<- se1[m]^2
          t1[m]<-summary(model1[[m]])$p.t[2]
          p1[m]<-summary(model1[[m]])$p.pv[2]
          rsq[m]<-summary(model1[[m]])$r.sq
          gcv[m]<-model1[[m]]$gcv
        }else{
          coef1[m]<-NA
          se1[m]<-NA
          var1[m]<- NA
          t1[m]<-NA
          p1[m]<-NA
          rsq[m]<-NA
          gcv[m]<-NA
        }
        if (!is.null(model2[[m]]) & !is.infinite(model2[[m]]$gcv.ubre)){
          coef2[m]<-summary(model2[[m]])$p.coeff[2]
          se2[m]<-summary(model2[[m]])$se[2]
          var2[m]<- se2[m]^2
          t2[m]<-summary(model2[[m]])$p.t[2]
          p2[m]<-summary(model2[[m]])$p.pv[2]
          rsq[m]<-summary(model2[[m]])$r.sq
          gcv[m]<-model2[[m]]$gcv
        }else{
          coef2[m]<-NA
          se2[m]<-NA
          var2[m]<- NA
          t2[m]<-NA
          p2[m]<-NA
          rsq[m]<-NA
          gcv[m]<-NA
        }
        code_city[m] <- list1[[m]]$code_city[1]
      }
      #print(paste0(Sys.time(), " model1 ", m))
      list3[[j]] <- data.table(outcome = j,strata = hw,code_city,coef1,se1,var1,t1,p1,coef2,se2,var2,t2,p2,rsq,gcv)
    }
    list2[[hw]] <- rbindlist(list3)
    rm(list3);gc()
  }
}


save(list2, file = "/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/list_pm25_hw_strata.Rdata")
data.pm25 <- rbindlist(list2)
fwrite(data.pm25, "/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/list_pm25_hw_strata.csv")

#meta
data.pm25 <- fread("/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/list_pm25_hw_strata.csv")
data.pm25[abs(coef1) > 0.1 | se1 > 0.1,]$var1 <- NA
data.pm25[abs(coef2) > 0.1 | se2 > 0.1,]$var2 <- NA

list3 <- list()
for (d in c("lag2","lag3")) {
  for (p in c("p90", "p92.5", "p95", "p97.5")) {
    list4 <- list()
    hw <- paste0("hw_mean_",d,"_",p)
    #hw <- "hw_mean_lag2_p90"
    for (i in c("total", "cvd", "res")){
      #i <- "total"
      a <- data.pm25[outcome == i, ]
      a <- a[strata == hw, ]
      print(paste0(Sys.time(), " meta ", i, " ", hw))
      indnomiss <- which(!is.na(a$var1))
      y1 <- a[indnomiss,coef1]
      S1 <- a[indnomiss,var1]
      meta1 <- mixmeta(y1, S1, random=~1|code_city, data=data.city[indnomiss,], 
                       control=list(showiter=F))
      beta1<-meta1$coefficients
      se1<-sqrt(meta1$vcov)
      betaoverall1<-round((exp(beta1*10)-1)*100,2)
      betaoverallinf1<-round((exp((beta1*10-1.96*se1*10))-1)*100,2)
      betaoverallsup1<-round((exp((beta1*10+1.96*se1*10))-1)*100,2)
      
      indnomiss <- which(!is.na(a$var2))
      y2 <- a[indnomiss,coef2]
      S2 <- a[indnomiss,var2]
      meta2 <- mixmeta(y2, S2, random=~1|code_city, data=data.city[indnomiss,], 
                       control=list(showiter=F))
      beta2<-meta2$coefficients
      se2<-sqrt(meta2$vcov)
      betaoverall2<-round((exp(beta2*10)-1)*100,2)
      betaoverallinf2<-round((exp((beta2*10-1.96*se2*10))-1)*100,2)
      betaoverallsup2<-round((exp((beta2*10+1.96*se2*10))-1)*100,2)
      
      z_strata <- (beta2 - beta1) / sqrt(se1^2 + se2^2)
      p_strata <- 2 * pnorm(-abs(z_strata))
      
      list4[[i]] <- data.table(outcome = i, strata = hw,
                               n_hw = paste0(betaoverall2, " (",betaoverallinf2,",",betaoverallsup2, ")"),
                               hw    = paste0(betaoverall1, " (",betaoverallinf1,",",betaoverallsup1, ")"))
      
      list4[[i]]$p <- p_strata 
    }
    list3[[hw]] <- rbindlist(list4)
    rm(list4);gc()
  }
}
pool.pm25 <- rbindlist(list3)
fwrite(pool.pm25, "/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/pool_pm25_hw_strata.csv")
rm(list3);gc()
    


