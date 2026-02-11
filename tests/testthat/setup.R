# tests/testthat/setup.R
# 这个文件会在所有测试前自动运行

# 生成测试数据
set.seed(123)
n <- 50
p <- 2
x <- seq(-5, 5, length.out = n)
X <- cbind(1, x, x^2)
theta_true <- c(-6, -3, 0.25)
e_sd <- 0.2

# 生成多个 Y 样本
Y.list <- list()
for (i in 1:10) {
  set.seed(1000 + i)
  epsilon <- rnorm(n, mean = 0, sd = e_sd)
  Y.list[[i]] <- X %*% theta_true + epsilon
}

# 这些对象在所有测试中都可用
