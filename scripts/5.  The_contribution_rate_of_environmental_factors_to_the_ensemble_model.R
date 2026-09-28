
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
    VarImport = 0,  
    modeling.id = "modtest"  
  )

# ============================================================================
# 集成模型：只保留一个最终集成模型，并计算环境因子贡献率
# ============================================================================

# 这里只构建 TSS 加权平均集成模型（weighted mean ensemble）。
# 不再进行当前/未来投影，也不绘制响应曲线。
Pafei_integrated_models <-
  BIOMOD_EnsembleModeling(
    modeling.output = Pafei_models,
    em.by = "all",
    eval.metric = "TSS",
    eval.metric.quality.threshold = 0.8,
    models.eval.meth = c("KAPPA", "TSS", "ROC"),
    prob.mean = FALSE,
    prob.cv = FALSE,
    prob.ci = FALSE,
    prob.median = FALSE,
    committee.averaging = FALSE,
    prob.mean.weight = TRUE,
    prob.mean.weight.decay = "proportional",
    VarImport = 10
  )

# 提取集成模型变量重要性。
# biomod2 3.5.1 返回数组：第1维为环境变量，其余维度为置换重复/集成模型。
Pafei_integrated_var_import <-
  get_variables_importance(Pafei_integrated_models)

# 对每个环境变量在所有置换重复上求平均，得到一套重要性数值。
var_import_mean <-
  apply(Pafei_integrated_var_import, 1, mean, na.rm = TRUE)

# 归一化为贡献率（%），使所有环境变量贡献率之和 = 100%。
var_contribution <- data.frame(
  Variable = names(var_import_mean),
  Contribution_percent = 100 * var_import_mean / sum(var_import_mean, na.rm = TRUE),
  row.names = NULL,
  check.names = FALSE
)

# 按贡献率从高到低排列，并保留两位小数。
var_contribution <- var_contribution[
  order(var_contribution$Contribution_percent, decreasing = TRUE),
]
var_contribution$Contribution_percent <-
  round(var_contribution$Contribution_percent, 2)

# 控制台只输出这一套最终贡献率。
print(var_contribution)
cat("贡献率合计 = ", sum(var_contribution$Contribution_percent), "%\n", sep = "")

# 导出 CSV。
write.csv(
  var_contribution,
  file = "Pafei_ensemble_variable_contribution_percent.csv",
  row.names = FALSE
)
