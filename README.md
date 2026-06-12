# Adaptive spatial blocking for scalable clustering inference

See the full paper on arXiv: https://arxiv.org/pdf/2606.12021

## Overview

This repository provides an R implementation of **B-KAMP**, a block-based framework for scalable spatial clustering inference in large-scale spatial point patterns.

The framework consists of two main components:

### `adaptive_spatial_blocking()`

This function identifies a set of rectangular blocks within a given observation window using a grid representation.

The procedure applies an adaptive spatial blocking algorithm to partition the observation window into blocks that satisfy point-count requirements and shape constraints.

This algorithm enables fast preprocessing for large-scale images and local spatial analysis.

### `bkamp_test()`

This function performs the proposed block-based spatial clustering test.

Instead of evaluating clustering evidence over the full observation window, the method detects local clustering signals within each block and aggregates the resulting evidence through a weighted test statistic.

By replacing a full-window pairwise computation with block-level computations, the proposed approach improves computational scalability while retaining sensitivity to local clustering signals.

## Arguments

### `adaptive_spatial_blocking()`

| Argument      | Description                                                              |
| ------------- | ------------------------------------------------------------------------ |
| `grid`        | A list of grid cells.                                                    |
| `n_min`       | Minimum point-count requirements for valid blocks.                       |
| `rho1`        | Maximum aspect ratio during Phase I.                                     |
| `rho2`        | Maximum aspect ratio during Phase II.                                    |
| `bivariate`   | Logical indicator for bivariate analysis.                                |
| `target_mark` | Mark corresponding to the point type of interest in univariate analysis. |
| `type1_mark`  | Mark corresponding to the first point type in bivariate analysis.        |
| `type2_mark`  | Mark corresponding to the second point type in bivariate analysis.       |

### `bkamp_test()`

| Argument    | Description                                                                                                               |
| ----------- | ------------------------------------------------------------------------------------------------------------------------- |
| `pp_list`   | A list of selected blocks returned by `adaptive_spatial_blocking()`.                                                      |
| `rratio`    | Radius ratio for K-function evaluation. The evaluation radius is defined as `rratio × min(block width, block height)`. |
| `bivariate` | Logical indicator for bivariate analysis.                                                                                 |

## Values

### `adaptive_spatial_blocking()`

| Value                 | Description                                                                                    |
| --------------------- | ---------------------------------------------------------------------------------------------- |
| `success_blocks`      | Identified blocks within the observation window that satisfy the block constraints. |
| `residual_grid_cells` | Remaining grid cells that are not assigned to any identified block.                            |
| `assigned_vector`     | Logical indicator showing whether each grid cell is assigned to a selected block.              |

### `bkamp_test()`

| Value       | Description                                                                  |
| ----------- | ---------------------------------------------------------------------------- |
| `ratio`     | Radius ratio used for K-function evaluation.                                 |
| `statistic` | Test statistic Z(r) obtained by aggregating block-level clustering evidence. |
| `p_value`   | Approximate p-value from the standard normal distribution.                   |
| `b`         | Number of identified blocks.                    |

## Example

The following example generates a marked spatial point pattern on the unit square with three point types, `"a"`, `"b"`, and `"c"`. Each point type is independently generated from a homogeneous Poisson point process with intensity 1000, 2000, and 3000, respectively.

```r
library(spatstat.geom)

set.seed(123)

W <- owin(c(0, 1), c(0, 1))

pp_a <- rpoispp(lambda = 1000, win = W)
marks(pp_a) <- rep("a", pp_a$n)

pp_b <- rpoispp(lambda = 2000, win = W)
marks(pp_b) <- rep("b", pp_b$n)

pp_c <- rpoispp(lambda = 3000, win = W)
marks(pp_c) <- rep("c", pp_c$n)

ppp <- superimpose(
  pp_a,
  pp_b,
  pp_c,
  W = W
)
```

### Univariate clustering test

The following code applies B-KAMP to test the univariate clustering of point type `"a"`.

```r
split_max <- floor(sqrt(ppp$n))

grid <- construct_adaptive_grid(
  ppp_obj   = ppp,
  split_max = split_max
)

target_mark <- "a"
m_target <- sum(marks(ppp) == target_mark)
m_bg     <- ppp$n - m_target

n_min <- c(
  sqrt(m_target),
  min(sqrt(m_target), sqrt(m_bg))
)

res_blocking <- adaptive_spatial_blocking(
  grid        = grid,
  n_min       = n_min,
  rho1        = 1,
  rho2        = Inf,
  bivariate   = FALSE,
  target_mark = target_mark
)

blocks <- res_blocking$success_blocks

res_uni <- bkamp_test(
  pp_list   = blocks,
  bivariate = FALSE,
  rratio    = seq(0.05, 0.2, 0.05)
)

res_uni
```

### Bivariate colocalization test

The following code applies B-KAMP to test the bivariate colocalization between point types `"a"` and `"b"`.

```r
split_max <- floor(sqrt(ppp$n))

grid <- construct_adaptive_grid(
  ppp_obj   = ppp,
  split_max = split_max
)

target1 <- "a"
target2 <- "b"

m_target1 <- sum(marks(ppp) == target1)
m_target2 <- sum(marks(ppp) == target2)
m_bg      <- ppp$n - (m_target1 + m_target2)

n_min <- c(
  sqrt(m_target1),
  sqrt(m_target2),
  min(sqrt(m_target1), sqrt(m_target2), sqrt(m_bg))
)

res_blocking <- adaptive_spatial_blocking(
  grid       = grid,
  n_min      = n_min,
  rho1       = 1,
  rho2       = Inf,
  bivariate  = TRUE,
  type1_mark = target1,
  type2_mark = target2
)

blocks <- res_blocking$success_blocks

res_biv <- bkamp_test(
  pp_list   = blocks,
  bivariate = TRUE,
  rratio    = seq(0.05, 0.2, 0.05)
)

res_biv
```

## Results

For the univariate clustering test of type `"a"`, B-KAMP returns the following result.

```r
# A tibble: 4 × 4
  ratio statistic p_value     b
  <dbl>     <dbl>   <dbl> <int>
1  0.05   -0.697    0.757    17
2  0.10   -0.675    0.750    17
3  0.15   -0.474    0.682    17
4  0.20   -0.0586   0.523    17
```

For the bivariate colocalization test between types `"a"` and `"b"`, B-KAMP returns the following result.

```r
# A tibble: 4 × 4
  ratio statistic p_value     b
  <dbl>     <dbl>   <dbl> <int>
1  0.05    -0.424  0.664     17
2  0.10    -0.202  0.580     17
3  0.15     0.892  0.186     17
4  0.20     1.78   0.0379    17
```
