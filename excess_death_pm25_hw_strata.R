### 272 excess death pm25_hw
### sunwr
### 2025/5/6

library(tsModel)
library(data.table)
library(mgcv)
library(splines)
library(dplyr)
library(zoo)
library(mixmeta)

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
list_hw <- lapply(list1, function(x) x[,heatwave:= yg.hw(temp, 2, 0.925)])


##pm25
a <- fread("/Users/sunweiran/Desktop/interaction_SWR/result/fenceng/pool_pm25_hw_strata.csv")
#272循环
beta0 <- 0.001646501 
se0 <- 0.0008059329
beta3 <- 0.01263518  # 热浪天RR的对数
se3 <- 0.002937537          # 热浪天RR对数标准误
beta4 <- 0.00147023  # 非热浪天RR的对数
se4 <- 0.0005463598        # 非热浪天RR对数标准误
RR_0 <- exp(beta0)
RR_h <- exp(beta3)   # 热浪天的RR值
RR_nh <- exp(beta4)  # 非热浪天的RR值
C0 <- 15         
delta_C_unit <- 10

results <- list()
for (i in 1:272) {
  df <- list_hw[[i]]
  print(i)
  df <- df %>% 
    select(date, total, heatwave, l02pm25) %>%
    na.omit()  

  # 计算每日超额死亡
  df <- df %>%
    mutate(
      delta_C = pmax(l02pm25 - C0, 0),
      
      RR_0 = RR_0^(delta_C / delta_C_unit), 
      AF_0 = (RR_0 - 1) / RR_0,
      excess_deaths_0 = AF_0 * total,
      
      
      RR = ifelse(heatwave == 1, RR_h, RR_nh),
      RR_i = RR^(delta_C / delta_C_unit),
      AF = (RR_i - 1) / RR_i,
      excess_deaths = AF * total
    )
  
  # 创建包含i的结果数据框
  result <- data.frame(
    id = i, # 添加标识符
    hw_nointer = sum(df$excess_deaths_0[df$heatwave == 1], na.rm = TRUE),  
    hw_inter = sum(df$excess_deaths[df$heatwave == 1], na.rm = TRUE),
    nohw_nointer = sum(df$excess_deaths_0[df$heatwave == 0], na.rm = TRUE),  
    nohw_inter = sum(df$excess_deaths[df$heatwave == 0], na.rm = TRUE),
    total_nointer = sum(df$excess_deaths_0, na.rm = TRUE),  
    total_inter = sum(df$excess_deaths, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # 将结果添加到结果列表
  results[[i]] <- result
}

# 将所有结果合并为一个数据框
final_results <- do.call(rbind, results)
result_1 <- data.frame(
  hw_nointer = sum(final_results$hw_nointer, na.rm = TRUE),  
  hw_inter = sum(final_results$hw_inter, na.rm = TRUE),
  nohw_nointer = sum(final_results$nohw_nointer, na.rm = TRUE),  
  nohw_inter = sum(final_results$nohw_inter, na.rm = TRUE),
  total_nointer = sum(final_results$total_nointer, na.rm = TRUE),  
  total_inter = sum(final_results$total_inter, na.rm = TRUE),
  stringsAsFactors = FALSE
)
fwrite(result_1, "/Users/sunweiran/Desktop/interaction_SWR/result/AFP/excess_death_pm25_hw.csv")


# 初始化存储所有模拟结果的列表
set.seed(22573)       # 设定随机种子保证可重复性
# 蒙特卡洛模拟次数
n_sim <- 10000
result_sim <- data.frame()
for (sim in 1:n_sim) {
  print(sim)
  # 生成随机RR值（对数正态分布）
  RR_0_sim <- exp(rnorm(1, mean = beta0, sd = se0))
  RR_h_sim <- exp(rnorm(1, mean = beta3, sd = se3))
  RR_nh_sim <- exp(rnorm(1, mean = beta4, sd = se4))
  
  # 存储当前模拟的272个城市结果
  sim_results <- list()
  
  for (i in 1:272) {
    df <- list_hw[[i]] %>% 
      select(date, total, heatwave, l02pm25) %>%
      na.omit()
    
    # 使用模拟的RR值计算超额死亡
    df <- df %>%
      mutate(
        delta_C = pmax(l02pm25 - C0, 0),
        
        # 无交互作用模型
        RR_0 = RR_0_sim^(delta_C / delta_C_unit),
        AF_0 = (RR_0 - 1) / RR_0,
        excess_deaths_0 = AF_0 * total,
        
        # 有交互作用模型（使用当前模拟的RR）
        RR = ifelse(heatwave == 1, RR_h_sim, RR_nh_sim),
        RR_i = RR^(delta_C / delta_C_unit),
        AF = (RR_i - 1) / RR_i,
        excess_deaths = AF * total
      )
    
    sim_results[[i]] <- data.frame(
      id = i,
      sim_id = sim,
      hw_nointer = sum(df$excess_deaths_0[df$heatwave == 1], na.rm = TRUE),
      hw_inter = sum(df$excess_deaths[df$heatwave == 1], na.rm = TRUE),
      nohw_nointer = sum(df$excess_deaths_0[df$heatwave == 0], na.rm = TRUE),
      nohw_inter = sum(df$excess_deaths[df$heatwave == 0], na.rm = TRUE),
      total_nointer = sum(df$excess_deaths_0, na.rm = TRUE),
      total_inter = sum(df$excess_deaths, na.rm = TRUE)
    )
  }
  # 合并当前模拟的所有城市结果
  final_sim_df <- do.call(rbind, sim_results)
  
  # 计算当前模拟的全国总和并添加到result_sim
  current_sim_result <- data.frame(
    sim = sim,
    hw_nointer = sum(final_sim_df$hw_nointer, na.rm = TRUE),  
    hw_inter = sum(final_sim_df$hw_inter, na.rm = TRUE),
    nohw_nointer = sum(final_sim_df$nohw_nointer, na.rm = TRUE),  
    nohw_inter = sum(final_sim_df$nohw_inter, na.rm = TRUE),
    total_nointer = sum(final_sim_df$total_nointer, na.rm = TRUE),  
    total_inter = sum(final_sim_df$total_inter, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  result_sim <- rbind(result_sim, current_sim_result)
}
conf_intervals <- sapply(result_sim[, -1], function(x) {
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
})

# 转置并格式化结果
result <- as.data.frame(t(conf_intervals))
colnames(result) <- c("95% CI Lower", "95% CI Upper")
result
fwrite(result_sim, "/Users/sunweiran/Desktop/interaction_SWR/result/AFP/excess_death_pm25_hw_10000.csv")
fwrite(result, "/Users/sunweiran/Desktop/interaction_SWR/result/AFP/results_excess_death_pm25_hw_10000.csv")

hw_p_value <- mean(result_sim$hw_inter < result_sim$hw_nointer)
nohw_p_value <- mean(result_sim$nohw_inter > result_sim$nohw_nointer)
total_p_value <- mean(result_sim$total_inter < result_sim$total_nointer)



p_value <- t.test(result_sim$hw_inter, result_sim$total_nointer, var.equal = TRUE)$p.value
p_value <- t.test(result_sim$nohw_inter, result_sim$nohw_nointer, var.equal = TRUE)$p.value
p_value <- t.test(result_sim$total_inter, result_sim$total_nointer, var.equal = TRUE)$p.value





