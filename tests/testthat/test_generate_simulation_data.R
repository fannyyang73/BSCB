test_that("generate_simulation_data returns correct structure", {
  sim <- generate_simulation_data(
    p            = 2,
    n            = 20,
    e_sd         = 0.2,
    theta_true   = c(-6, -3, 0.25),
    replication  = 2,
    design_index = 2
  )

  expect_type(sim, "list")
  expect_named(sim, c("X", "Y.list", "optimal_x", "optimal_weights"))
})

test_that("X has correct dimensions", {
  sim <- generate_simulation_data(
    p = 2, n = 20, e_sd = 0.2,
    theta_true = c(-6, -3, 0.25),
    replication = 2
  )
  expect_equal(nrow(sim$X), 20)
  expect_equal(ncol(sim$X), 3)   # p + 1 = 3
})

test_that("Y.list has correct length and each Y has length n", {
  sim <- generate_simulation_data(
    p = 2, n = 20, e_sd = 0.2,
    theta_true = c(-6, -3, 0.25),
    replication = 2
  )
  expect_length(sim$Y.list, 2)
  expect_true(all(sapply(sim$Y.list, length) == 20))
})

test_that("results are reproducible with same batch_index and same seed", {
  sim1 <- generate_simulation_data(
    p = 2, n = 20, e_sd = 0.2,
    theta_true = c(-6, -3, 0.25),
    replication = 2, batch_index = 1, seed=42
  )
  sim2 <- generate_simulation_data(
    p = 2, n = 20, e_sd = 0.2,
    theta_true = c(-6, -3, 0.25),
    replication = 2, batch_index = 1, seed=42
  )
  expect_identical(sim1$Y.list, sim2$Y.list)
})

test_that("different batch_index produces different Y", {
  sim1 <- generate_simulation_data(
    p = 2, n = 20, e_sd = 0.2,
    theta_true = c(-6, -3, 0.25),
    replication = 2, batch_index = 1
  )
  sim2 <- generate_simulation_data(
    p = 2, n = 20, e_sd = 0.2,
    theta_true = c(-6, -3, 0.25),
    replication = 2, batch_index = 2
  )
  expect_false(identical(sim1$Y.list[[1]], sim2$Y.list[[1]]))
})

test_that("design_index = 1 works and weights sum to n", {
  sim <- generate_simulation_data(
    p = 2, n = 20, e_sd = 0.2,
    theta_true = c(-6, -3, 0.25),
    replication = 2, design_index = 1
  )
  expect_equal(sum(sim$optimal_weights), 20)
})

test_that("invalid inputs are rejected", {
  expect_error(
    generate_simulation_data(p = 4, n = 20, e_sd = 0.2,
                             theta_true = c(1, 2, 3, 4, 5), replication = 2),
    "p must be 1, 2, or 3"
  )
  expect_error(
    generate_simulation_data(p = 2, n = 20, e_sd = 0.2,
                             theta_true = c(1, 2), replication = 2),
    "length of theta_true must be p\\+1"
  )
})
