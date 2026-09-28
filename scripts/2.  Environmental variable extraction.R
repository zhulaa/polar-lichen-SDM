setwd("F:/Part 3_Prediction of lichen distribution/tripolar_maxent/An/env")

# 1. 加载必要的 R 包
# 安装必要的包
install.packages(c("raster", "sf", "terra"))
# 加载包
library(raster)
library(terra)
library(sf)


# 2. 批量导入环境数据

env_data_path <- "J:/Climatedata/Climatedata_TIFF_30s" # 指定存放环境数据的文件夹路径
env_files <- list.files(env_data_path, pattern = "\\.tif$", full.names = TRUE) # 获取所有 .tif 文件的完整路径
env_rasters <- lapply(env_files, function(file) {
  raster <- rast(file)  # 加载栅格
  names(raster) <- tools::file_path_sans_ext(basename(file))  # 使用文件名作为图层名称
  return(raster)
})

lapply(env_rasters, names)

# 批量导入所有 .tif 文件并堆叠为多层栅格对象
env_stack <- rast(env_rasters)  # 使用raster包的stack函数
print(env_stack) # 查看栅格数据的基本信息
names(env_stack) # 查看栅格图层名称
plot(env_stack[[1]], main = names(env_stack)[1]) # 查看某一层数据


# 3. 加载研究区矢量文件
# 由于rgdal包依赖于4.3.*以下的R版本，因此，用sf替代
install.packages("rnaturalearth")
install.packages("remotes")  # 确保 remotes 包已安装
remotes::install_github("ropensci/rnaturalearthhires") # 高分辨率扩展
library(rnaturalearth)
library(rnaturalearthhires)

# 获取南极洲矢量数据
antarctica <- ne_countries(scale = "large", continent = "Antarctica", returnclass = "sf")

# 将 sf 对象转换为 SpatVector（terra 使用的矢量格式）
study_area_vect <- vect(antarctica)

# 查看和绘制
plot(study_area_vect)


# 4. 裁剪环境数据至研究区
# 裁剪包括两步：裁剪范围和掩膜处理。
env_crop <- crop(env_stack, study_area_vect) # 裁剪栅格数据到研究区范围
env_masked <- mask(env_crop, study_area_vect) # 使用研究区掩膜裁剪后的栅格

# 检查裁剪结果
plot(env_masked[[5]], main = "Cropped Environmental Data")


# 5. 检查和标准化数据
# 检查分辨率和对齐：确保所有栅格层具有相同的分辨率和范围。如果发现不一致，可以使用 terra::resample 进行对齐：
template_layer <- env_masked[[1]] # 使用第一层作为模板对齐其他图层
env_aligned <- resample(env_masked, template_layer)

# 保存裁剪后的栅格数据
writeRaster(env_aligned, "env_aligned_30s.tif", 
            filetype = "GTiff",  # 保存为 GeoTIFF 格式
            overwrite = TRUE)  # 如果已存在文件，允许覆盖

#------------------------------------顺利运行到此-------------------------------

# 去除多重共线性 = Pearson相关性系数进行筛选，去除高度相关的变量，避免多重共线性对建模的影响。
# 方法1：计算皮尔森相关性矩阵，可视化相关性矩阵，根据阈值筛选变量(可用作方法2的验证)
# 将环境数据转为数据框
env_values_pearson <- as.data.frame(terra::extract(env_stack, terra::xyFromCell(env_stack, seq_len(ncell(env_stack)))))
cor_matrix <- cor(env_values_pearson, use = "complete.obs", method = "pearson") # 计算相关性矩阵

print(cor_matrix) # 查看相关性矩阵
write.csv(cor_matrix, "correlation_matrix.csv", row.names = TRUE) # 导出相关性矩阵为 CSV 文件

#  可视化相关性矩阵：可以通过热图更直观地查看变量之间的相关性：
library(corrplot)
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.cex = 0.7)
#手动筛选：通过查阅相关性矩阵或热图，选择每对高度相关的变量中保留一个进行建模。

# 方法2： 使用包 usdm 提供的 VIF（方差膨胀因子） 方法，更系统地移除多余变量。usdm 会基于多重共线性自动筛选出较优的变量集合。
install.packages("usdm")
library(usdm)

# 将环境变量转换为数据框检查VIF
env_values <- as.data.frame(values(env_aligned))
colnames(env_values) <- make.names(names(env_aligned),unique = TRUE) # 确保列名唯一且无重复
env_values_VIF <- vifstep(env_values, th = 10) # 设置 VIF 阈值（如 10）
reduced_variables <- exclude(env_values, env_values_VIF)  # 排除具有高共线性的变量

print(env_values_VIF)


#----------------------------占内存过大，运行以下内容---------------------------

## 内存过大时！
print(env_aligned)  # 检查栅格数据的大小和分辨率
terra::ncell(env_aligned)  # 查看总像元数
terra::memory.size()  # 检查当前 R 进程内存占用（仅 Windows）


# 分块处理栅格数据 (内存过大)
library(terra)
output_file <- "env_values_partial.csv"
if (file.exists(output_file)) file.remove(output_file)
write.table(data.frame(x = NA, y = NA, layer1 = NA), output_file, 
            sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)

# 分块处理并逐块保存
for (block in blocks) {
  block_coords <- terra::xyFromCell(env_aligned, block)   # 坐标
  block_values <- terra::extract(env_aligned, block)      # 环境值
  
  # 合并坐标与环境值
  block_data <- cbind(block_coords, block_values)
  
  # 逐块追加保存
  write.table(block_data, output_file, sep = ",", 
              row.names = FALSE, col.names = FALSE, append = TRUE)
}

# 检查文件
print(paste("Data saved to:", output_file))


block_size <- 1e6 # 定义块大小（例如，每块提取 10,000,000 个像元的值）
n_cells <- terra::ncell(env_aligned)
blocks <- split(1:n_cells, ceiling(seq_along(1:n_cells) / block_size)) # 分块索引
env_values <- NULL # 初始化结果

# 遍历每块并提取数据
for (block in blocks) {
  # 获取当前块的坐标和栅格值
  block_coords <- terra::xyFromCell(env_aligned, block)  # 获取坐标
  block_values <- terra::extract(env_aligned, block_coords)  # 提取环境值
  
  # 合并当前块的数据
  env_values <- rbind(env_values, cbind(block_coords, block_values))
  
  # 可选：保存中间结果到磁盘，减少内存占用
  write.csv(env_values, "env_values_partial.csv", row.names = FALSE, append = TRUE)
}

# 查看结果
head(env_values)








