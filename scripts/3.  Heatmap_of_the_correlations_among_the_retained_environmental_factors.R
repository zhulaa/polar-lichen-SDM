


Environmental_data <- list.files("J:/BIOMOD2/biomod2/Antarctica/EnvironmentData/current", pattern = ".asc$", full.names = TRUE) %>%
  stack()
Environmental_data <- list.files("J:/BIOMOD2/biomod2/Arctic/EnvironmentData/current", pattern = ".asc$", full.names = TRUE) %>%
  stack()
Environmental_data


Environmental_data_df <- na.omit(as.data.frame(Environmental_data))
head(Environmental_data_df)

#相关性热图
Environmental_data_df %>%
  as.matrix() %>%
  quickcor(cor.test = TRUE) +
  geom_square(data = get_data(type = "upper", show.diag = FALSE)) +
  geom_number(aes(num = r),  size = 5, family = "Times New Roman",
              data = get_data(type = "lower", show.diag = FALSE)) +
  geom_abline(slope = -1, linewidth = 0.8, intercept = 24) +
  scale_fill_gradientn(colours = c("#00AFBB", "white",  "#FC4E07"))+
  
  theme_minimal() + # 使用简洁的主题
  theme(
    
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    # 修改坐标轴标题
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10), family = "Times New Roman"), # X轴标题
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(t = 50), family = "Times New Roman"), # Y轴标题
    # 修改坐标轴刻度标签
    axis.text.x = element_text(size = 14, color = "black", hjust = 0.5, family = "Times New Roman"),  # X轴刻度标签
    axis.text.y = element_text(size = 14, color = "black", family = "Times New Roman"), # Y轴刻度标签
    # 修改图例标题和标签（图例字体为新罗马）
    legend.title = element_text(size = 12, family = "Times New Roman"),  # 图例标题字体（新罗马）
    legend.text = element_text(size = 10, family = "Times New Roman")  # 图例标签字体（新罗马）
    
  ) +
  
  labs(
    x = "Environmental variable",
    y = "Environmental variable"
  )+
  
  theme(plot.margin = margin(t = 20, b = 20))  # 增加上、下的边距



