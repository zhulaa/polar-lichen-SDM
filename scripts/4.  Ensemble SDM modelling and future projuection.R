
knitr::opts_chunk$set(
  echo = FALSE,
  fig.align = "center",
  message = FALSE,
  warning = FALSE,
  comment = "#>",
  dpi = 1200
)




setwd("F:/BIOMOD2/biomod2/Antarctica/SpeciesData1")

# 安装biomod2包
if (!require(remotes)) install.packages("remotes")
if (!require(biomod2)) remotes::install_version("biomod2", version = "3.5.1")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(raster)) install.packages("raster")
if (!require(tmap)) install.packages("tmap")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(sf)) install.packages("sf")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(raster)) install.packages("raster")
if (!require(terra)) install.packages("terra")
if (!require(tidyterra)) install.packages("tidyterra")
if (!require(factoextra)) install.packages("factoextra")
if (!require(ggfortify)) install.packages("ggfortify")
if (!require(ade4)) install.packages("ade4")
if (!require(ggcor)) remotes::install_github("Github-Yilei/ggcor")
if (!require(sjmisc)) install.packages("sjmisc")
if (!require(flextable)) install.packages("flextable")
if (!require(xlsx)) install.packages("xlsx")

library(cli)
setwd("F:/BIOMOD2/biomod2/Antarctica/SpeciesData1")
#导入分布点和环境变量数据
spe_example <- read.csv("spe_pale.csv") 
head(spe_example)



Environmental_data <- list.files("F:/BIOMOD2/biomod2/Antarctica/EnvironmentData/current", pattern = ".asc$", full.names = TRUE) %>%
  stack()
Environmental_data



#创建伪不存在点

Pafei_data <-
  BIOMOD_FormatingData(
    resp.var = spe_example["Presences.Absences"],
    resp.xy = spe_example[, c("longitude", "latitude")],
    expl.var = Environmental_data,  
    resp.name = "Presences.Absences", 
    PA.nb.rep = 2,              
    PA.nb.absences = 1000,      
    PA.strategy = "random"
  )
Pafei_data    
#结果显示的有以下信息

plot(Pafei_data)



#设置模型参数

Model_parameter <-
  BIOMOD_ModelingOptions(
    GLM = list(
      type = "polynomial",
      interaction.level = 0,
      myFormula = NULL,
      test = "AIC",
      family = binomial(link = "logit"),
      control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE),
      mustart = 0.5
    ),
    #GLM广义线性模型
    GBM = list(
      distribution = "bernoulli",
      n.trees = 2000,
      interaction.depth = 7,
      n.minobsinnode = 5,
      shrinkage = 0.001,
      bag.fraction = 0.5,
      train.fraction = 1,
      keep.data = FALSE,
      verbose = FALSE,
      perf.method = "cv",
      n.cores = 8,
      cv.folds = 3
    ),
    #GBM梯度提升模型
    GAM = list(
      algo = "GAM_mgcv",
      type = "s_smoother",
      k = -1,
      interaction.level = 0,
      myFormula = NULL,
      family = binomial(link = "logit"),
      method = "GCV.Cp",
      select = FALSE,
      knots = NULL,
      paraPen = NULL,
      control = list(
        nthreads = 16,
        irls.reg = 0,
        epsilon = 1e-07,
        maxit = 200,
        trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15, rank.tol = 1.49011611938477e-08,
        nlm = list(ndigit = 7, gradtol = 1e-06, stepmax = 2, steptol = 1e-04, iterlim = 200, check.analyticals = 0),
        optim = list(factr = 1e+07), newton = list(conv.tol = 1e-06, maxNstep = 5, maxSstep = 2, maxHalf = 30, use.svd = 0),
        outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, efs.lspmax = 15,
        efs.tol = 0.1, keepData = FALSE, scale.est = "fletcher", edge.correct = FALSE
      ),
      optimizer = c("outer", "newton")
    ),
    #GAM广义相加模型
    CTA = list(
      method = "class",
      parms = "default",
      cost = NULL,
      control = list(
        xval = 5,
        minbucket = 5,
        minsplit = 5,
        cp = 0.001,
        maxdepth = 25
      )
    ),
    #CTA模型
    ANN = list(
      NbCV = 5,
      size = NULL,
      decay = NULL,
      rang = 0.1,
      maxit = 200
    ),
    SRE = list(quant = 0.025),
    FDA = list(
      method = "mars",
      add_args = NULL
    ),
    #ANN人工神经网络模型
    MARS = list(
      type = "simple",
      interaction.level = 0,
      myFormula = NULL,
      nk = NULL,
      penalty = 2,
      thresh = 0.001,
      nprune = NULL,
      pmethod = "backward"
    ),
    #MARS多元回归适应样条
    RF = list(
      do.classif = TRUE,
      ntree = 500,
      mtry = "default",
      nodesize = 5,
      maxnodes = NULL
    ),
    #RF随机森林模型
    MAXENT.Phillips = list(
      path_to_maxent.jar = "F:/BIOMOD2/biomod2/Antarctica/maxent.jar",    #maxent主程序的存放路径要修改
      memory_allocated = 6144,
      background_data_dir = "default",
      maximumbackground = "default",
      maximumiterations = 1000,
      visible = FALSE,
      linear = TRUE,
      quadratic = TRUE,
      product = TRUE,
      threshold = TRUE,
      hinge = TRUE,
      lq2lqptthreshold = 80,
      l2lqthreshold = 10,
      hingethreshold = 15,
      beta_threshold = -1,
      beta_categorical = -1,
      beta_lqp = -1,
      beta_hinge = -1,
      betamultiplier = 1,
      defaultprevalence = 0.5
      #MAXENT最大熵模型
    )
  )


# 模型运行（一般文章会选择5-6种模型运行）YAFEI:DATAFORM:http://lucky-boy.ysepan.com/
# ANN: 人工神经网络 Artificial neural net work; 
#CART: 分类回归树分析 Classification and regression tree analysis; 
#FDA: 柔性判别分析 Flexible discriminant analisis; 
#GAM: 广义相加模型Generalized additive model; 
#GBM: 助推法 Generalized boosting model; 
#GLM: 广义线性模型 Generalized linear models; 
#MARS: 多元自适应回归样条模型 Multiple adaptive regression splines; 
#MAXENT: 最大熵模型 Maximum entropy models; 
#RF: 随机森林 Radom Forest; 
#SRE: 表面分布区室模型 Surface range envelope



Pafei_models<-
  BIOMOD_Modeling(
    data = Pafei_data,    
    models = c(
      "GLM", "GBM", "CTA", "ANN", "SRE", "FDA", "RF", "GAM", "MARS", "MAXENT.Phillips"
    ),             # "GLM", "GBM", "GAM", "CTA", "ANN", "SRE", "FDA", "MARS", "RF", "MAXENT.Phillips"
    models.options = Model_parameter,  
    NbRunEval = 2,    
    DataSplit = 75,   
    VarImport = 4,  
    modeling.id = "modtest"  
  )




#结果评价分析 

Pafei_models_scores <- get_evaluations(Pafei_models)
dim(Pafei_models_scores)
dimnames(Pafei_models_scores)

Pafei_models_scores


#模型评估指标转换表格并导出表格

Pafei_models_scores<-Pafei_models_scores %>%
  as.data.frame.array() %>%
  sjmisc::rotate_df(rn = "data")

data <- Pafei_models_scores

file_path <- "F:/BIOMOD2/biomod2/Antarctica/SpeciesData1/Pafei_models_scores（模型评估指标）.csv"

write.csv(data, file = file_path, row.names = FALSE)


#-------------------------模型评价结果-------------------------------
# 定义指标和模型
metrics <- c("TSS", "AUC", "KAPPA")  # 评估指标
models <- unique(Pafei_models_scores$Model)  # 模型名称

# 初始化结果存储
summary_stats <- data.frame()

# 循环计算每个模型和指标的统计量
for (metric in metrics) {
  for (model in models) {
    # 筛选当前模型和指标的值
    metric_values <- Pafei_models_scores[Pafei_models_scores$Model == model, metric]
    
    # 计算平均值、标准差、变异系数
    mean_value <- mean(metric_values, na.rm = TRUE)
    sd_value <- sd(metric_values, na.rm = TRUE)
    cv_value <- (sd_value / mean_value) * 100  # 变异系数(%)
    
    # 添加到结果表
    summary_stats <- rbind(summary_stats, data.frame(
      Model = model,
      Metric = metric,
      Mean = mean_value,
      SD = sd_value,
      CV = cv_value
    ))
  }
}

# 查看结果
print(summary_stats)

# 保存结果到CSV
write.csv(summary_stats, "Pafei_models_evaluation_summary.csv", row.names = FALSE)


#----------------------------------






#下面的6张图无法绘制！！！！！！！！！！！！
#绘制"ROC"-"TSS"指标图
pdf("ROC-TSS指标图1.pdf")   
models_scores_graph(
  Pafei_models,
  by = "models",
  metrics = c("ROC", "TSS"),
  xlim = c(0, 1),  
  ylim = c(0, 1)   
)



pdf("ROC-TSS指标图2.pdf")   
models_scores_graph(
  Pafei_models,
  by = "cv_run",
  metrics = c("ROC", "TSS"),
  xlim = c(0, 1),   
  ylim = c(0, 1)   
)



pdf("ROC-TSS指标图3.pdf")   
models_scores_graph(
  Pafei_models,
  by = "data_set",
  metrics = c("ROC", "TSS"),
  xlim = c(0, 1),   
  ylim = c(0, 1)   
)


#绘制"ROC"-"KAPPA"指标图


pdf("ROC-KAPPA指标图1.pdf")   
models_scores_graph(
  Pafei_models,
  by = "models",
  metrics = c("ROC", "KAPPA"),
  xlim = c(0, 1),   
  ylim = c(0, 1)   
)



pdf("ROC-KAPPA指标图2.pdf")   
models_scores_graph(
  Pafei_models,
  by = "cv_run",
  metrics = c("ROC", "KAPPA"),
  xlim = c(0, 1),   
  ylim = c(0, 1)   
)



pdf("ROC-KAPPA指标图3.pdf")   
models_scores_graph(
  Pafei_models,
  by = "data_set",
  metrics = c("ROC", "KAPPA"),
  xlim = c(0, 1),   
  ylim = c(0, 1)   
)



###变量重要性(权重)


Pafei_models_var_import <- get_variables_importance(Pafei_models)
Pafei_models_var_import

#将上面的结果结果转表格导出，如果模型重复了两次不要运行这一步，直接运行下一个不然报错


Pafei_models_var_import1<-Pafei_models_var_import%>%
  as.data.frame.array() %>%
  sjmisc::rotate_df(rn = "data")

# 读取Pafei_models_scores数据集
data <- Pafei_models_var_import1

# 设置导出文件路径和文件名
file_path <- "J:/BIOMOD2/biomod2/Antarctica/SpeciesData1/Pafei_models_var_import（变量重要性权重）.csv"

# 使用write.csv()函数将数据导出为CSV文件
write.csv(data, file = file_path, row.names = FALSE)




###计算均值（模型未重复的话默认跳过该步骤）


apply(Pafei_models_var_import, c(1, 2), mean) %>%    #设置了两次以上创建伪不存在点，即存在PA1和PA2时运行改脚本，求两次的平均值
  as.data.frame() %>%
  rownames_to_column("vars") %>%
  flextable::regulartable()   #将结果转换为表格




###响应曲线


Pafei_glm <- BIOMOD_LoadModels(Pafei_models, models="GLM")    
Pafei_glm


#响应曲线图中有多条线

pdf("GLM模型响应曲线2D图.pdf")   
myRespPlot2D <- 
  response.plot2(
    models = Pafei_glm,
    Data = get_formal_data(Pafei_models, 'expl.var'),
    show.variables = get_formal_data(Pafei_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var.metric = 'median',
    col = c("red", "blue", "green", "purple", "orange", "cyan"),
    legend = TRUE,
    data_species = get_formal_data(Pafei_models, 'resp.var')
  )
dim(myRespPlot2D)   #响应曲线绘制的原始数据
dimnames(myRespPlot2D)



###3D响应变量图

pdf("GLM模型响应曲线3D图.pdf")
myRespPlot3D <- 
  response.plot2(
    models = Pafei_glm[1],
    Data = get_formal_data(Pafei_models, 'expl.var'), 
    show.variables = get_formal_data(Pafei_models, 'expl.var.names'),
    do.bivariate = TRUE,
    fixed.var.metric = 'median',
    data_species = get_formal_data(Pafei_models, 'resp.var'),
    display_title = FALSE
  )
dim(myRespPlot3D) 
dimnames(myRespPlot3D)


### 集成模型（重要）



Pafei_integrated_models <-
  BIOMOD_EnsembleModeling(
    modeling.output = Pafei_models,
    em.by = "all",            
    eval.metric = "TSS",     
    eval.metric.quality.threshold = 0.8,     
    models.eval.meth = c("KAPPA","TSS", "ROC"),  
    prob.mean = TRUE, 
    prob.cv = TRUE,   
    committee.averaging = TRUE,      #相对多数投票法
    prob.mean.weight = TRUE, #权重概率法
    VarImport = 0
  )


###集成模型评估


Pafei_integrated_models_scores <- get_evaluations(Pafei_integrated_models)
Pafei_integrated_models_scores



###集成模型预测


Pafei_models_proj_current <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data,
    proj.name = "current",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_models_proj_current


Pafei_integrated_models_proj_current <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_current,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_current



#---------------新增!!!!!!!!!!!!!  2025.1.22

# 获取已训练模型的名称
model_names <- BIOMOD_LoadModels(Pafei_models, models = "GLM")  # 只加载GLM模型的名称

# 获取所有气候因子变量名称
climate_variables <- get_formal_data(Pafei_models, 'expl.var.names')

# 定义颜色向量（根据运行结果）
run_colors <- c("red", "blue", "green", "purple", "orange", "cyan")

# 遍历每个气候因子，绘制响应曲线并保存为PDF文件
for (var in climate_variables) {
  # 定义输出PDF文件名
  pdf_filename <- paste0(var, "_响应曲线2D图.pdf")
  
  # 打开PDF设备
  pdf(pdf_filename)
  
  # 绘制响应曲线
  response.plot2(
    models = model_names,  # 使用模型名称字符向量
    Data = get_formal_data(Pafei_models, 'expl.var'),  # 环境变量数据
    show.variables = var,  # 当前气候因子
    do.bivariate = FALSE,  # 单变量响应曲线
    fixed.var.metric = 'median',  # 其他变量固定为中位数
    col = run_colors,  # 设置运行结果的颜色
    legend = TRUE,  # 是否显示图例
    data_species = get_formal_data(Pafei_models, 'resp.var')  # 物种分布数据
  )
  
  # 关闭PDF设备
  dev.off()
}




###新增2025.3.31





#----------------------------已运行至此，没有问题-------------------------------
#--------------------------如报错，请重新打开并运行-----------------------------




#---
### 9.预测未来(未来时期目录下的环境数据个数必须与当前一致，并且一一对应)

### 加载2050年数据


#绘制2050环境数据图
pdf("2050环境数据.pdf")   
Environmental_data_2050_RCP45 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2050SSP126", pattern = ".asc$", full.names = TRUE) %>%
  stack()
Environmental_data
tm_shape(Environmental_data_2050_SSP126) +
  tm_raster(title = "value") +
  tm_facets(ncol = 4)


### 1. 对2050SSP126预测
Environmental_data_2050_SSP126 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2050SSP126", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2050_SSP126 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2050_SSP126,
    proj.name = "2050_SSP126",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2050_SSP126 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2050_SSP126,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2050_SSP126


### 预测效果
pdf("2050年SSP126集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2050_SSP126, str.grep = "EMca|EMwmean")


### 2. 对2050SSP245预测
Environmental_data_2050_SSP245 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2050SSP245", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2050_SSP245 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2050_SSP245,
    proj.name = "2050_SSP245",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2050_SSP245 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2050_SSP245,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2050_SSP245


### 预测效果
pdf("2050年SSP245集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2050_SSP245, str.grep = "EMca|EMwmean")


### 3. 对2050SSP370预测
Environmental_data_2050_SSP370 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2050SSP370", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2050_SSP370 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2050_SSP370,
    proj.name = "2050_SSP370",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2050_SSP370 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2050_SSP370,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2050_SSP370


### 预测效果
pdf("2050年SSP370集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2050_SSP370, str.grep = "EMca|EMwmean")



### 4. 对2050SSP585预测
Environmental_data_2050_SSP585 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2050SSP585", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2050_SSP585 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2050_SSP585,
    proj.name = "2050_SSP585",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2050_SSP585 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2050_SSP585,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2050_SSP585


### 预测效果
pdf("2050年SSP585集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2050_SSP585, str.grep = "EMca|EMwmean")



-----------------------------------
  ### 加载2090数据
  
  
  #绘制2090环境数据图
  pdf("2090环境数据.pdf")   
Environmental_data_2070_RCP45 <- list.files("F:/YF2024/EnvironmentData/2070RCP45", pattern = ".asc$", full.names = TRUE) %>%
  stack()
Environmental_data
tm_shape(Environmental_data_2050_RCP45) +
  tm_raster(title = "value") +
  tm_facets(ncol = 4
  )


### 对2090预测

### 5. 对2090SSP126预测
Environmental_data_2090_SSP126 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2090SSP126", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2090_SSP126 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2090_SSP126,
    proj.name = "2090_SSP126",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2090_SSP126 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2090_SSP126,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2090_SSP126


### 预测效果
pdf("2090年SSP126集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2090_SSP126, str.grep = "EMca|EMwmean")


### 6. 对2090SSP245预测
Environmental_data_2090_SSP245 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2090SSP245", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2090_SSP245 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2090_SSP245,
    proj.name = "2090_SSP245",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2090_SSP245 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2090_SSP245,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2090_SSP245


### 预测效果
pdf("2090年SSP245集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2090_SSP245, str.grep = "EMca|EMwmean")



### 7. 对2090SSP370预测
Environmental_data_2090_SSP370 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2090SSP370", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2090_SSP370 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2090_SSP370,
    proj.name = "2090_SSP370",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2090_SSP370 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2090_SSP370,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2090_SSP370


### 预测效果
pdf("2090年SSP370集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2090_SSP370, str.grep = "EMca|EMwmean")




### 8. 对2090SSP585预测
Environmental_data_2090_SSP585 <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/future/SSP/2090SSP585", pattern = ".asc$", full.names = TRUE) %>%
  stack()

Pafei_models_proj_2090_SSP585 <-
  BIOMOD_Projection(
    modeling.output = Pafei_models,
    new.env = Environmental_data_2090_SSP585,
    proj.name = "2090_SSP585",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
#下面用集成模型预测未来
Pafei_integrated_models_proj_2090_SSP585 <-
  BIOMOD_EnsembleForecasting(
    EM.output = Pafei_integrated_models,
    projection.output = Pafei_models_proj_2090_SSP585,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
Pafei_integrated_models_proj_2090_SSP585


### 预测效果
pdf("2090年SSP585集成模型EMca-EMwm预测结果.pdf")   
plot(Pafei_integrated_models_proj_2090_SSP585, str.grep = "EMca|EMwmean")




### 10.范围变化



Pafei_bin_proj_current <-
  raster::stack(
    c(
      ca = "Presences.Absences/proj_current/individual_projections/Presences.Absences_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
      wm = "Presences.Absences/proj_current/individual_projections/Presences.Absences_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
    )
  )


Pafei_bin_proj_2050_RCP45 <-
  raster::stack(
    c(
      ca = "Presences.Absences/proj_2050_RCP45/individual_projections/Presences.Absences_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
      wm = "Presences.Absences/proj_2050_RCP45/individual_projections/Presences.Absences_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
    )
  )

Pafei_bin_proj_2070_RCP45 <-
  raster::stack(
    c(
      ca = "Presences.Absences/proj_2070_RCP45/individual_projections/Presences.Absences_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
      wm = "Presences.Absences/proj_2070_RCP45/individual_projections/Presences.Absences_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
    )
  )


### 当前 -> 2050


SRC_current_2050_RCP45 <-
  BIOMOD_RangeSize(
    Pafei_bin_proj_current,
    Pafei_bin_proj_2050_RCP45
  )
SRC_current_2050_RCP45$Compt.By.Models %>% 
  as.data.frame.array() %>% 
  rownames_to_column("type") %>% 
  regulartable()


### 当前 -> 2070


SRC_current_2070_RCP45 <-
  BIOMOD_RangeSize(
    Pafei_bin_proj_current,
    Pafei_bin_proj_2070_RCP45
  )
SRC_current_2070_RCP45$Compt.By.Models %>% 
  as.data.frame.array() %>% 
  rownames_to_column("type") %>% 
  regulartable()


### 适宜区变化制图


pdf("当前到未来的适宜区变化-EMca-EMw变化图.pdf")   
Pafei_src_map <-
  raster::stack(
    SRC_current_2050_RCP45$Diff.By.Pixel,
    SRC_current_2070_RCP45$Diff.By.Pixel
  )
names(Pafei_src_map) <- c("ca cur-2050", "wm cur-2050", "ca cur-2070", "wm cur-2070")
my.at <- seq(-2.5, 1.5, 1)
myColorkey <-
  list(
    at = my.at, 
    labels =
      list(
        labels = c("lost", "pres", "abs", "gain"), 
        at = my.at[-1] - 0.5 
      )
  )
rasterVis::levelplot(
  Pafei_src_map,
  main = "Species Pafei range change",
  colorkey = myColorkey,
  col.regions = c("#f03b20", "#99d8c9", "#f0f0f0", "#2ca25f"),
  layout = c(2, 2)
)


#---
