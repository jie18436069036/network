#清除所有变量
rm(list = ls())
#################### 创建必要的函数######################
# 计算向量中零的数量
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)}
# 计算向量中负值的平均值
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
  }
# 计算向量中正值的平均值
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
  }
############## 工作流程选项########################
## 设置分类持续性阈值（最小存在比例），用于保留分析中的分类
pers.cutoff <- 0.10
## 设置虚拟模型的迭代次数（建议>=200）
iter <- 200
## 选择分类洗牌（tax.shuffle = T）或行洗牌（tax.shuffle= F）
tax.shuffle <- T
## 选择是否使用自定义相关矩阵
# 注意：您的相关性表必须与丰度表有相同数量的分类。丰度表中不应有空的（全零）分类向量。
# 即使您输入了自定义的相关性表，持续性阈值也会被应用
use.custom.cors <- F
########### 计算内聚力#######################
## 读取数据集（每行是一个样本）
b <-
  read.csv("RH_3h.csv", header = T, row.names = 1)
# 若使用自定义相关矩阵，请读取并检查维度
if(use.custom.cors == T) {
  custom.cor.mat <- read.csv("RH_3h.csv", header = T,
                             row.names = 1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  print(dim(b)[2] == dim(custom.cor.mat)[2])
  }
# 重新格式化数据，去除空样本和分类
c <- as.matrix(b)
c <- c[rowSums(c) > 0,
       colSums(c) > 0]
# 保存原始样本的总个体数
rowsums.orig <- rowSums(c)
# 根据持续性阈值确定分类零的数量阈值
zero.cutoff <- ceiling(pers.cutoff
                       * dim(c)[1])
# 删除低于持续性阈值的分类
d <- c[ , apply(c, 2, zero) <
          (dim(c)[1]-zero.cutoff) ]
# 删除没有个体的样本
d <- d[rowSums(d) > 0, ]
# 如果使用自定义相关矩阵，更新相关矩阵
if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) <
                                         (dim(c)[1]-zero.cutoff), apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
  }
# 创建相对丰度矩阵
rel.d <- d / rowsums.orig
# 查看保留的社区比例
hist(rowSums(rel.d))
# 计算观察到的相关矩阵
cor.mat.true <- cor(rel.d)
# 保存每个分类的中位数相关性
med.tax.cors <- vector()
# 计算虚拟模型的预期相关性# 如果使用自定义相关矩阵，则跳过虚拟模型
if(use.custom.cors == F) {
  if(tax.shuffle) {
    for(which.taxon in 1:dim(rel.d)[2]){
# 保存每次置换的相关性
      perm.cor.vec.mat <- vector()
for(i in 1:iter){
  # 创建与 rel.d 维度相同的空矩阵
  perm.rel.d <- matrix(numeric(0),
                       dim(rel.d)[1], dim(rel.d)[2])
  rownames(perm.rel.d) <-
    rownames(rel.d)
  colnames(perm.rel.d) <-
    colnames(rel.d)
  # 对每个分类进行洗牌
  for(j in 1:dim(rel.d)[2]){
    perm.rel.d[, j ] <- sample(rel.d[
      ,j ])
    }
  # 保持焦点列不变
  perm.rel.d[, which.taxon] <- rel.d[
    , which.taxon]
  # 计算置换矩阵的相关性矩阵
  cor.mat.null <- cor(perm.rel.d)
  # 保存焦点分类的相关性
  perm.cor.vec.mat <-
    cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
}
      # 保存中位数相关性
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1,
                                                median))
if(which.taxon %% 20 == 0){print(which.taxon)}
      }
    } else {
      for(which.taxon in 1:dim(rel.d)[2]){
  # 保存每次置换的相关性
        perm.cor.vec.mat <- vector()
  for(i in 1:iter){
    # 复制矩阵以进行丰度洗牌
    perm.rel.d <- rel.d
    # 对每个样本进行洗牌
    for(j in 1:dim(rel.d)[1]){
      which.replace <- which(rel.d[j, ]
                             > 0 )
      which.replace.nonfocal <-
        which.replace[!(which.replace %in% which.taxon)]
    # 置换分类向量
      perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal])
      }
    # 计算置换矩阵的相关性矩阵
    cor.mat.null <- cor(perm.rel.d)
    # 保存焦点分类的相关性
    perm.cor.vec.mat <-
      cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
  }
        # 保存中位数相关性
        med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1,
                                                  median))
  if(which.taxon %% 20 == 0){print(which.taxon)}
        }
      }
  }
# 计算观察到的减去预期的相关性
if(use.custom.cors == T) {
  obs.exp.cors.mat <- custom.cor.mat.sub
  } else {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
    }
diag(obs.exp.cors.mat) <- 0
######## 生成连通度和内聚度向量 # 计算连通度（正负相关性平均值）
connectedness.pos <-
  apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <-
  apply(obs.exp.cors.mat, 2, neg.mean)
# 计算内聚度（相对丰度与连通度相乘）
cohesion.pos <- rel.d %*%
  connectedness.pos
cohesion.neg <- rel.d %*%
  connectedness.neg
######## 输出结果 
output <- list(connectedness.neg,
               connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness","Positive Connectedness", "Negative Cohesion","Positive Cohesion")
print(output)
#####输出Connectedness数据矩阵############
# 提取要合并的矩阵或向量
connectedness_data <- list(
  "Negative Connectedness" = output[["Negative
                                     Connectedness"]],
  "Positive Connectedness" = output[["Positive
                                     Connectedness"]]
  )
# 转换为数据框，并确保行数相同
connectedness_df <-
  as.data.frame(do.call(cbind, connectedness_data))
# 为每列设置名称
colnames(connectedness_df) <-
  names(connectedness_data)
# 写入 CSV 文件
write.csv(connectedness_df, file ="connectedness.csv", row.names = TRUE)
###############输出Cohesion数据###########
# 提取要合并的矩阵或向量
cohesion_data <- list(
  "Negative Cohesion" = output[["Negative Cohesion"]],
  "Positive Cohesion" = output[["Positive Cohesion"]])
# 转换为数据框，并确保行数相同
cohesion_df <-
  as.data.frame(do.call(cbind, cohesion_data))
# 为每列设置名称
colnames(cohesion_df) <-
names(cohesion_data)
# 写入 CSV 文件
write.csv(cohesion_df, file ="cohesion_RH_3h.csv", row.names = TRUE)
