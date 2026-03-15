# 加载必要的库
library(igraph)

# 读取CSV文件
file_path <- "227edge_RH_3h.csv"  # 修改为你的CSV文件路径
df <- read.csv(file_path, stringsAsFactors = FALSE)

# 创建图
G <- graph_from_data_frame(d = df, directed = FALSE)

# 计算自然连通度的函数
natural_connectivity <- function(G) {
  # 获取图的邻接矩阵
  A <- as_adjacency_matrix(G, sparse = FALSE)
  
  # 计算邻接矩阵的特征值
  eigenvalues <- eigen(A)$values
  
  # 确保特征值是非负的
  eigenvalues <- Re(eigenvalues)  # 取实部
  eigenvalues <- eigenvalues[eigenvalues > 0]  # 去掉负值和零值
  
  # 计算自然连通度
  n <- vcount(G)  # 图的节点数
  if (n == 0) {
    return(0.0)
  }
  phi <- sum(eigenvalues) / n
  return(phi)
}

# 定义要移除的节点比例
proportions <- c(0.1, 0.2, 0.3, 0.4, 0.5)  # 10%, 20%, 30%, 40%, 50%
results <- data.frame(ProportionRemoved = proportions, NaturalConnectivity = numeric(length(proportions)))

# 计算移除指定比例节点后的自然连通度
for (i in seq_along(proportions)) {
  p <- proportions[i]
  # 随机移除指定比例的节点
  nodes_to_remove <- sample(V(G), size = round(p * vcount(G)))
  G_reduced <- delete_vertices(G, nodes_to_remove)
  results$NaturalConnectivity[i] <- natural_connectivity(G_reduced)
}

# 输出结果
print(results)

# 保存结果到CSV文件
output_file_path <- "natural_connectivity_results_RH_3h.csv"  # 输出文件路径
write.csv(results, file = output_file_path, row.names = FALSE)
