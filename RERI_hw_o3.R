library(tsModel)
library(data.table)
library(mgcv)
library(splines)
library(dplyr)
library(zoo)
library(mixmeta)
library(MASS)

load("/Users/sunweiran/Desktop/interaction_SWR/raw_data/272city_0216.Rdata")
#na.fill(ts,"extend")

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

#####fixed_O3_cut_offs#####
list2 <- list()
for (a in c("l02o3")) {
  list3 <- list()
  list.cat1 <- lapply(list_hw, function(x) x[, cat.o3 := as.numeric(cut(get(a), c(0,100,999)))])
    for (b in c("mean")){
      for (d in c("lag2","lag3")){
        for (p in c("p90", "p92.5", "p95", "p97.5")) {
          list.cat2 <- lapply(list.cat1, function(x) x[, cat.hw := get(paste0("hw_",b,"_",d,"_",p))+1])
          list.cat2 <- lapply(list.cat2, function(x) x[, cat := ifelse((is.na(cat.o3) | is.na(cat.hw)) == T, NA, paste0(cat.o3, cat.hw))])
          list4 <- list()
          for (j in c("total", "cvd", "res")) {
            fml <- get(j)~cat+ns(l02rh,3)+as.factor(dow)+ns(date,df=round(7*length(date)/365.25))
            model1 <- list()
            for(i in 1:272) {
              print(paste0(Sys.time(), " ", j, " ", i, " ",b," ",d," ",p, " ", a))
              dt <- list.cat2[[i]]
              tryCatch({
                model1[[i]] <- gam(fml, data = dt, family=quasipoisson,na.action=na.exclude)
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
            
            coef3<-vector()
            se3<-vector()
            var3<-vector()
            t3<-vector()
            p3<-vector()
            
            rsq<-vector()
            gcv<-vector()
            
            cov12<-vector()
            cov13<-vector()
            cov23<-vector()
            
            code_city<-vector()
            
            nmod <- length(model1)
            if (nmod>0){
              for (m in 1:nmod){
                if (!is.null(model1[[m]] )){
                  if (!is.infinite(model1[[m]]$gcv.ubre)){
                    coef1[m]<-summary(model1[[m]])$p.coeff[2]
                    se1[m]<-summary(model1[[m]])$se[2]
                    var1[m]<- se1[m]^2
                    t1[m]<-summary(model1[[m]])$p.t[2]
                    p1[m]<-summary(model1[[m]])$p.pv[2]
                    coef2[m]<-summary(model1[[m]])$p.coeff[3]
                    se2[m]<-summary(model1[[m]])$se[3]
                    var2[m]<- se2[m]^2
                    t2[m]<-summary(model1[[m]])$p.t[3]
                    p2[m]<-summary(model1[[m]])$p.pv[3]
                    rsq[m]<-summary(model1[[m]])$r.sq
                    gcv[m]<-model1[[m]]$gcv
                    if(names(summary(model1[[m]])$p.coeff)[4]=="cat22"){
                      coef3[m]<-summary(model1[[m]])$p.coeff[4]
                      se3[m]<-summary(model1[[m]])$se[4]
                      var3[m]<- se3[m]^2
                      t3[m]<-summary(model1[[m]])$p.t[4]
                      p3[m]<-summary(model1[[m]])$p.pv[4] 
                      # 提取三组别系数(cat12/cat21/cat22)的 3x3 城市内协方差(Bayesian后验)
                      Vm <- tryCatch(vcov(model1[[m]], freq = FALSE), error = function(e) NULL)
                      if (!is.null(Vm) && nrow(Vm) >= 4) {
                        V3 <- Vm[2:4, 2:4]
                        cov12[m] <- V3[1, 2]; cov13[m] <- V3[1, 3]; cov23[m] <- V3[2, 3]
                      }
                    }
                    else{
                      coef3[m]<-NA
                      se3[m]<-NA
                      var3[m]<- NA
                      t3[m]<-NA
                      p3[m]<-NA 
                    }
                  } else{
                    coef1[m]<-NA
                    se1[m]<-NA
                    var1[m]<- NA
                    t1[m]<-NA
                    p1[m]<-NA
                    coef2[m]<-NA
                    se2[m]<-NA
                    var2[m]<- NA
                    t2[m]<-NA
                    p2[m]<-NA
                    coef3[m]<-NA
                    se3[m]<-NA
                    var3[m]<- NA
                    t3[m]<-NA
                    p3[m]<-NA
                    rsq[m]<-NA
                    gcv[m]<-NA
                  }
                }
                else{
                  coef1[m]<-NA
                  se1[m]<-NA
                  var1[m]<- NA
                  t1[m]<-NA
                  p1[m]<-NA
                  coef2[m]<-NA
                  se2[m]<-NA
                  var2[m]<- NA
                  t2[m]<-NA
                  p2[m]<-NA
                  coef3[m]<-NA
                  se3[m]<-NA
                  var3[m]<- NA
                  t3[m]<-NA
                  p3[m]<-NA
                  rsq[m]<-NA
                  gcv[m]<-NA
                }
                #print(paste0(Sys.time(), " model1 ", m))
                code_city[m] <- list1[[m]]$code_city[1]
              }
              
            }
            list4[[j]] <- data.table(outcome = j, cut_o3=100,cut_hw = p, day_hw = d, tm = b,lag = a,
                                     code_city,coef1,se1,var1,t1,p1,coef2,se2,var2,t2,p2,coef3,se3,var3,t3,p3,cov12,cov13,cov23,rsq,gcv)
          }
        list3[[paste0(b,"_",d,"_",p)]] <- rbindlist(list4)
        rm(list4);gc()
        }
      }
  }
  list2[[a]] <- rbindlist(list3)
  rm(list3);gc()
}
data <- rbindlist(list2)
save(list2, file = "/Users/sunweiran/Desktop/interaction_SWR/result/o3/list_o3_combine_hw_reri.Rdata")
fwrite(data, "/Users/sunweiran/Desktop/interaction_SWR/result/o3/data_o3_combine_hw_reri.csv")

#meta
data <- fread("/Users/sunweiran/Desktop/interaction_SWR/result/o3/data_o3_combine_hw_reri.csv")
data[abs(coef1) > 10 | se1 > 10, c("var1", "cov12", "cov13") := .(NA, NA, NA)]
data[abs(coef2) > 10 | se2 > 10, c("var2", "cov12", "cov23") := .(NA, NA, NA)]
data[abs(coef3) > 10 | se3 > 10, c("var3", "cov13", "cov23") := .(NA, NA, NA)]

list6 <- list()
for (j in c("total", "cvd", "res")) {
  for (l in c("l02o3")) {
      for (b in c("mean")) {
        for (d in c("lag2","lag3")) {
          for (p in c("p90", "p92.5", "p95", "p97.5")) {
            a <- data[outcome==j & cut_hw==p & day_hw==d & lag==l & tm==b,]
            print(paste0(Sys.time(), " ", j, " ", b," ", d," ", p, " ", l))
            #RR12
            indnomiss1 <- which(!is.na(a$var1))
            y1 <- a[indnomiss1,coef1]
            S1 <- a[indnomiss1,var1]
            meta1 <- mixmeta(y1, S1, random=~1|code_city, data=data.city[indnomiss1,], 
                             control=list(showiter=F))
            beta1<-meta1$coefficients
            se1<-sqrt(meta1$vcov)
            betaoverall1<-exp(beta1)
            betaoverallinf1<-exp(beta1-1.96*se1)
            betaoverallsup1<-exp(beta1+1.96*se1)
            #RR21
            indnomiss2 <- which(!is.na(a$var2))
            y2 <- a[indnomiss2,coef2]
            S2 <- a[indnomiss2,var2]
            meta2 <- mixmeta(y2, S2, random=~1|code_city, data=data.city[indnomiss2,], 
                             control=list(showiter=F))
            beta2<-meta2$coefficients
            se2<-sqrt(meta2$vcov)
            betaoverall2<-exp(beta2)
            betaoverallinf2<-exp(beta2-1.96*se2)
            betaoverallsup2<-exp(beta2+1.96*se2)
            ##RR22
            indnomiss3 <- which(!is.na(a$var3))
            y3 <- a[indnomiss3,coef3]
            S3 <- a[indnomiss3,var3]
            meta3 <- mixmeta(y3, S3, random=~1|code_city, data=data.city[indnomiss3,], 
                             control=list(showiter=F))
            beta3<-meta3$coefficients
            se3<-sqrt(meta3$vcov)
            betaoverall3<-exp(beta3)
            betaoverallinf3<-exp(beta3-1.96*se3)
            betaoverallsup3<-exp(beta3+1.96*se3)
            
            ## CI of RERI&SI：协方差修正（多元随机效应 meta，三系数联合合并）
            RERI_pt <- exp(beta3) - exp(beta2) - exp(beta1) + 1
            SI_pt   <- (exp(beta3) - 1) / ((exp(beta2) - 1) + (exp(beta1) - 1))
            ind_joint <- which(!is.na(a$var1) & !is.na(a$var2) & !is.na(a$var3) &
                                 !is.na(a$cov12) & !is.na(a$cov13) & !is.na(a$cov23))
            n_city_joint <- length(ind_joint)
            ci_RERI <- c(NA, NA, NA); ci_SI <- c(NA, NA, NA)
            if (n_city_joint > 5) {
              yJ <- as.matrix(a[ind_joint, .(coef1, coef2, coef3)])
              SJ <- lapply(seq_along(ind_joint), function(k) {
                r <- a[ind_joint[k], .(var1, cov12, cov13, var2, cov23, var3)]
                matrix(c(r$var1, r$cov12, r$cov13,
                         r$cov12, r$var2, r$cov23,
                         r$cov13, r$cov23, r$var3), 3, 3, byrow = TRUE)
              })
              metaJ <- tryCatch(mixmeta(yJ, SJ, random = ~1 | code_city, data = a[ind_joint, ],
                                        control = list(showiter = FALSE)), error = function(e) NULL)
              if (!is.null(metaJ)) {
                muJ <- as.numeric(coef(metaJ))
                VJ  <- vcov(metaJ)
                set.seed(20240118)
                draws <- tryCatch(mvrnorm(10000, mu = muJ, Sigma = VJ), error = function(e) NULL)
                if (!is.null(draws)) {
                  RERI_d <- exp(draws[, 3]) - exp(draws[, 2]) - exp(draws[, 1]) + 1
                  SI_d   <- (exp(draws[, 3]) - 1) / ((exp(draws[, 2]) - 1) + (exp(draws[, 1]) - 1))
                  ci_RERI <- quantile(RERI_d, c(0.025, 0.5, 0.975))
                  ci_SI   <- quantile(SI_d, c(0.025, 0.5, 0.975))
                  RERI_pt <- exp(muJ[3]) - exp(muJ[2]) - exp(muJ[1]) + 1
                  SI_pt   <- (exp(muJ[3]) - 1) / ((exp(muJ[2]) - 1) + (exp(muJ[1]) - 1))
                }
              }
            }
            
            list6[[paste0(j, " ", b," ", d," ", p, " ", l)]]<-data.frame(outcome = j,
                                                                         cut_pm25=100,
                                                                         cut_hw = p,
                                                                         day_hw = d,
                                                                         tm = b,
                                                                         lag = l,
                                   RR12 = paste0(round((betaoverall1-1)*100,2), " (",round((betaoverallinf1-1)*100,2),",",round((betaoverallsup1-1)*100,2), ")"),
                                   RR21 = paste0(round((betaoverall2-1)*100,2), " (",round((betaoverallinf2-1)*100,2),",",round((betaoverallsup2-1)*100,2), ")"),
                                   RR22 = paste0(round((betaoverall3-1)*100,2), " (",round((betaoverallinf3-1)*100,2),",",round((betaoverallsup3-1)*100,2), ")"),
                                   RERI = RERI_pt,
                                   si = SI_pt,
                                   ci_RERI = paste0(round(ci_RERI[2]*100,2), " (",round(ci_RERI[1]*100,2),",",round(ci_RERI[3]*100,2), ")"),
                                   ci_SI = paste0(round(ci_SI[2],3), " (",round(ci_SI[1],3),",",round(ci_SI[3],3), ")"), n_city_joint = n_city_joint
            )
        }
      }
    }
  }
}
pool.data <- rbindlist(list6)
fwrite(pool.data, "/Users/sunweiran/Desktop/interaction_SWR/result/o3/pool_o3_combine_hw_reri.csv")
