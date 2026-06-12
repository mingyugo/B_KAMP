required_pkgs <- c("spatstat.geom", "spatstat.explore", "dplyr", "purrr", "tibble")
invisible(lapply(required_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but is not installed.", pkg), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))


# ppp_obj: A marked point pattern of class `ppp`.
# split_max: Maximum allowed number of points in each grid cell.
# max_grid_dim: Maximum grid dimension to try.
construct_adaptive_grid <- function(ppp_obj, split_max, max_grid_dim = 100) {
  
  if (ppp_obj$n == 0) return(list())
  if (split_max < 1) stop("split_max must be at least 1.")
  
  win <- Window(ppp_obj)
  
  for (N in 2:max_grid_dim) {
    
    counts <- quadratcount(ppp_obj, nx = N, ny = N)
    
    max_cell_count <- max(as.vector(counts))
    
    if (max_cell_count <= split_max) {
      
      grid <- list()
      x_breaks <- seq(win$xrange[1], win$xrange[2], length.out = N + 1)
      y_breaks <- seq(win$yrange[1], win$yrange[2], length.out = N + 1)
      
      for (i in 1:N) {
        for (j in 1:N) {
          cell_win <- owin(
            xrange = c(x_breaks[j], x_breaks[j+1]),
            yrange = c(y_breaks[i], y_breaks[i+1])
          )
          grid[[length(grid) + 1]] <- ppp_obj[cell_win]
        }
      }
      return(grid)
    }
  }
  
  warning(sprintf("Could not find a grid up to %dx%d that satisfies the condition. Returning NULL.",
                  max_grid_dim, max_grid_dim))
  return(NULL)
}




# grid: A list of `ppp` objects returned by `construct_adaptive_grid()`.
# n_min: Minimum point-count constraints for valid blocks.
# rho1: Maximum aspect-ratio threshold for Phase I.
# rho2: Maximum aspect-ratio threshold for Phase II.
# bivariate: Logical; whether to use bivariate block constraints.
# target_mark: Mark value for the target point type in univariate analysis.
# type1_mark: Mark value for the first target point type in bivariate analysis.
# type2_mark: Mark value for the second target point type in bivariate analysis.
adaptive_spatial_blocking <- function(grid, n_min, rho1 = 1, rho2 = Inf,
                                      bivariate = FALSE,
                                      target_mark = "target",
                                      type1_mark = "target1",
                                      type2_mark = "target2",
                                      verbose = TRUE) {
  
  grid <- lapply(
    grid,
    standardize_bkamp_marks,
    bivariate = bivariate,
    target_mark = target_mark,
    type1_mark = type1_mark,
    type2_mark = type2_mark
  )
  
  ref_window <- spatstat.geom::Window(grid[[1]])
  ref_wd <- diff(ref_window$xrange)
  ref_hi <- diff(ref_window$yrange)
  base_tol <- sqrt(.Machine$double.eps) * 10
  
  assigned <- sapply(
    grid,
    check_block_constraints,
    n_min = n_min,
    bivariate = bivariate
  )
  
  success_blocks <- list()
  for (ind in which(assigned)) {
    success_blocks[[length(success_blocks) + 1]] <- list(ppp = grid[[ind]], ids = ind)
  }
  
  phase1_results <- run_shape_finding_phase(
    max_aspect_ratio = rho1,
    grid = grid,
    n_min = n_min,
    success_blocks = success_blocks,
    assigned = assigned,
    ref_wd = ref_wd,
    ref_hi = ref_hi,
    base_tol = base_tol,
    bivariate = bivariate,
    verbose = verbose
  )
  success_blocks <- phase1_results$success_blocks
  assigned <- phase1_results$assigned
  
  phase2_results <- run_shape_finding_phase(
    max_aspect_ratio = rho2,
    grid = grid,
    n_min = n_min,
    success_blocks = success_blocks,
    assigned = assigned,
    ref_wd = ref_wd,
    ref_hi = ref_hi,
    base_tol = base_tol,
    bivariate = bivariate,
    verbose = verbose
  )
  success_blocks <- phase2_results$success_blocks
  assigned <- phase2_results$assigned
  
  repeat {
    merge_made_in_this_pass <- FALSE
    leftover_ids <- which(!assigned)
    
    if (length(leftover_ids) == 0) break
    
    for (orphan_idx in leftover_ids) {
      local_candidate_merges <- list()
      
      for (success_idx in seq_along(success_blocks)) {
        orphan_win <- spatstat.geom::Window(grid[[orphan_idx]])
        success_win <- spatstat.geom::Window(success_blocks[[success_idx]]$ppp)
        
        if (length(get_adjacent_directions(orphan_win, success_win, tol = base_tol)) == 0) {
          next
        }
        
        cand_bbox <- spatstat.geom::boundingbox(orphan_win, success_win)
        
        ids_needed_to_fill <- which(sapply(grid, function(b) {
          w <- spatstat.geom::Window(b)
          (w$xrange[1] >= cand_bbox$xrange[1] - base_tol) &&
            (w$xrange[2] <= cand_bbox$xrange[2] + base_tol) &&
            (w$yrange[1] >= cand_bbox$yrange[1] - base_tol) &&
            (w$yrange[2] <= cand_bbox$yrange[2] + base_tol)
        }))
        
        current_success_ids <- success_blocks[[success_idx]]$ids
        valid_ids_for_merge <- c(current_success_ids, leftover_ids)
        
        if (!all(ids_needed_to_fill %in% valid_ids_for_merge)) next
        
        sub_blocks <- grid[ids_needed_to_fill]
        x_range <- range(unlist(lapply(sub_blocks, function(x) x$window$xrange)))
        y_range <- range(unlist(lapply(sub_blocks, function(x) x$window$yrange)))
        combined_window <- spatstat.geom::owin(xrange = x_range, yrange = y_range)
        
        merged_for_check <- do.call(
          spatstat.geom::superimpose,
          c(sub_blocks, list(W = combined_window))
        )
        
        if (!are_windows_equivalent(cand_bbox, spatstat.geom::Window(merged_for_check), tol = base_tol)) {
          next
        }
        
        new_width <- diff(cand_bbox$xrange)
        new_height <- diff(cand_bbox$yrange)
        
        score_ratio <- min(new_width, new_height) / max(new_width, new_height)
        score_points <- merged_for_check$n
        
        local_candidate_merges[[length(local_candidate_merges) + 1]] <- list(
          score_ratio = score_ratio,
          score_points = score_points,
          success_idx_to_update = success_idx,
          final_ids = ids_needed_to_fill,
          final_ppp = merged_for_check
        )
      }
      
      if (length(local_candidate_merges) > 0) {
        best_score_ratio <- -1
        best_score_points <- Inf
        best_candidate_idx <- -1
        
        for (i in seq_along(local_candidate_merges)) {
          cand <- local_candidate_merges[[i]]
          
          if (cand$score_ratio > best_score_ratio + base_tol) {
            best_score_ratio <- cand$score_ratio
            best_score_points <- cand$score_points
            best_candidate_idx <- i
          } else if (abs(cand$score_ratio - best_score_ratio) < base_tol) {
            if (cand$score_points < best_score_points) {
              best_score_points <- cand$score_points
              best_candidate_idx <- i
            }
          }
        }
        
        best_merge <- local_candidate_merges[[best_candidate_idx]]
        s_idx <- best_merge$success_idx_to_update
        success_blocks[[s_idx]]$ppp <- best_merge$final_ppp
        success_blocks[[s_idx]]$ids <- best_merge$final_ids
        assigned[best_merge$final_ids] <- TRUE
        
        merge_made_in_this_pass <- TRUE
        break
      }
    }
    
    if (!merge_made_in_this_pass) break
  }
  
  residual_grid_cells <- grid[which(!assigned)]
  
  if (verbose) {
    cat(sprintf(
      "\n--- Algorithm finished. Total %d success blocks found. ---\n",
      length(success_blocks)
    ))
  }
  
  return(list(
    success_blocks = lapply(success_blocks, function(x) x$ppp),
    residual_grid_cells = residual_grid_cells,
    assigned_vector = assigned
  ))
}

### helper functions for adaptive_spatial_blocking() ###
standardize_bkamp_marks <- function(pp, bivariate = FALSE,
                                    target_mark = "target",
                                    type1_mark = "target1",
                                    type2_mark = "target2") {
  mk <- as.character(spatstat.geom::marks(pp))
  
  if (bivariate) {
    mk[!mk %in% c(type1_mark, type2_mark)] <- "background"
    mk[mk == type1_mark] <- "target1"
    mk[mk == type2_mark] <- "target2"
  } else {
    mk[mk != target_mark] <- "background"
    mk[mk == target_mark] <- "target"
  }
  spatstat.geom::marks(pp) <- mk
  pp
}


check_block_constraints <- function(pp, n_min, bivariate = FALSE) {
  if (!bivariate) {
    target <- sum(pp$marks == "target")
    bg <- sum(pp$marks == "background")
    return((target >= n_min[1]) & (bg >= n_min[2]))
  }
  
  target1 <- sum(pp$marks == "target1")
  target2 <- sum(pp$marks == "target2")
  bg <- sum(pp$marks == "background")
  
  return((target1 >= n_min[1]) &
           (target2 >= n_min[2]) &
           (bg >= n_min[3]))
}


are_windows_equivalent <- function(win_target, win_actual, tol) {
  unfilled_area <- spatstat.geom::area.owin(
    spatstat.geom::setminus.owin(win_target, win_actual)
  )
  unfilled_area <= tol
}


get_adjacent_directions <- function(win_a, win_b, tol) {
  dirs <- c()
  
  overlap_y <- (min(win_a$yrange[2], win_b$yrange[2]) -
                  max(win_a$yrange[1], win_b$yrange[1])) > -tol
  overlap_x <- (min(win_a$xrange[2], win_b$xrange[2]) -
                  max(win_a$xrange[1], win_b$xrange[1])) > -tol
  
  if (abs(win_a$xrange[1] - win_b$xrange[2]) < tol && overlap_y) {
    dirs <- c(dirs, "left")
  }
  if (abs(win_a$xrange[2] - win_b$xrange[1]) < tol && overlap_y) {
    dirs <- c(dirs, "right")
  }
  if (abs(win_a$yrange[1] - win_b$yrange[2]) < tol && overlap_x) {
    dirs <- c(dirs, "down")
  }
  if (abs(win_a$yrange[2] - win_b$yrange[1]) < tol && overlap_x) {
    dirs <- c(dirs, "up")
  }
  
  unique(dirs)
}


find_best_rectangle <- function(grid_mat, max_aspect_ratio, ppp_obj, n_min,
                                x_coords, y_coords, ref_wd, ref_hi,
                                bivariate = FALSE) {
  if (sum(grid_mat) == 0) return(NULL)
  
  nr <- nrow(grid_mat)
  nc <- ncol(grid_mat)
  heights <- rep(0, nc)
  
  min_total_points <- Inf
  best_rect <- NULL
  
  for (i in seq_len(nr)) {
    heights <- ifelse(grid_mat[i, ] == 1, heights + 1, 0)
    stack <- list()
    
    for (j in seq_len(nc + 1)) {
      h <- if (j <= nc) heights[j] else 0
      
      while (length(stack) > 0 && heights[stack[[length(stack)]]$idx] >= h) {
        popped <- stack[[length(stack)]]
        stack[[length(stack)]] <- NULL
        
        height <- popped$h
        width <- if (length(stack) == 0) j - 1 else j - stack[[length(stack)]]$idx - 1
        
        if (height == 0 || width == 0) next
        
        aspect_ratio <- max(width * ref_wd, height * ref_hi) / min(width * ref_wd, height * ref_hi)
        if (aspect_ratio > max_aspect_ratio) next
        
        start_col <- if (length(stack) == 0) 1 else stack[[length(stack)]]$idx + 1
        x_start <- x_coords[start_col]
        x_end <- x_coords[start_col + width]
        r_start_in_grid <- i - height + 1
        y_start <- y_coords[nr - (r_start_in_grid + height) + 2]
        y_end <- y_coords[nr - r_start_in_grid + 2]
        
        cand_win <- spatstat.geom::owin(xrange = c(x_start, x_end), yrange = c(y_start, y_end))
        points_in_cand <- ppp_obj[cand_win]
        
        if (!check_block_constraints(points_in_cand, n_min = n_min, bivariate = bivariate)) next
        
        if (points_in_cand$n < min_total_points) {
          min_total_points <- points_in_cand$n
          best_rect <- c(row = r_start_in_grid, col = start_col, h = height, w = width)
        }
      }
      
      stack[[length(stack) + 1]] <- list(idx = j, h = h)
    }
  }
  
  return(best_rect)
}


run_shape_finding_phase <- function(max_aspect_ratio, grid, n_min,
                                    success_blocks, assigned,
                                    ref_wd, ref_hi, base_tol,
                                    bivariate = FALSE,
                                    verbose = TRUE) {
  global_win <- do.call(spatstat.geom::boundingbox, lapply(grid, spatstat.geom::Window))
  resolution_col <- min(sapply(grid, function(b) diff(spatstat.geom::Window(b)$xrange)))
  resolution_row <- min(sapply(grid, function(b) diff(spatstat.geom::Window(b)$yrange)))
  
  ncols <- round(diff(global_win$xrange) / resolution_col)
  nrows <- round(diff(global_win$yrange) / resolution_row)
  x_coords <- seq(from = global_win$xrange[1], to = global_win$xrange[2], length.out = ncols + 1)
  y_coords <- seq(from = global_win$yrange[1], to = global_win$yrange[2], length.out = nrows + 1)
  
  grid_mat <- matrix(1, nrow = nrows, ncol = ncols)
  s_wins <- lapply(success_blocks, function(list) spatstat.geom::Window(list$ppp))
  
  for (win in s_wins) {
    c_start <- floor((win$xrange[1] - global_win$xrange[1]) / resolution_col) + 1
    c_end <- floor((win$xrange[2] - global_win$xrange[1] - base_tol) / resolution_col) + 1
    r_start <- floor((win$yrange[1] - global_win$yrange[1]) / resolution_row) + 1
    r_end <- floor((win$yrange[2] - global_win$yrange[1] - base_tol) / resolution_row) + 1
    
    c_start <- max(1, c_start)
    c_end <- min(ncols, c_end)
    r_start <- max(1, r_start)
    r_end <- min(nrows, r_end)
    
    if (c_start <= c_end && r_start <= r_end) {
      grid_mat[(nrows - r_end + 1):(nrows - r_start + 1), c_start:c_end] <- 0
    }
  }
  
  unassigned_ids_current_phase <- which(!assigned)
  if (length(unassigned_ids_current_phase) == 0) {
    return(list(success_blocks = success_blocks, assigned = assigned))
  }
  
  sub_blocks <- grid[unassigned_ids_current_phase]
  x_range <- range(unlist(lapply(sub_blocks, function(x) x$window$xrange)))
  y_range <- range(unlist(lapply(sub_blocks, function(x) x$window$yrange)))
  combined_window <- spatstat.geom::owin(xrange = x_range, yrange = y_range)
  
  unassigned_ppp <- do.call(
    spatstat.geom::superimpose,
    c(sub_blocks, list(W = combined_window))
  )
  
  found_blocks <- list()
  
  repeat {
    candidate_info <- find_best_rectangle(
      grid_mat = grid_mat,
      max_aspect_ratio = max_aspect_ratio,
      ppp_obj = unassigned_ppp,
      n_min = n_min,
      x_coords = x_coords,
      y_coords = y_coords,
      ref_wd = ref_wd,
      ref_hi = ref_hi,
      bivariate = bivariate
    )
    
    if (is.null(candidate_info)) break
    
    r_final <- candidate_info[1]
    c_final <- candidate_info[2]
    h_final <- candidate_info[3]
    w_final <- candidate_info[4]
    
    x_start <- x_coords[c_final]
    x_end <- x_coords[c_final + w_final]
    y_start <- y_coords[nrows - (r_final + h_final - 1) + 1]
    y_end <- y_coords[nrows - r_final + 2]
    
    new_found_win <- spatstat.geom::owin(xrange = c(x_start, x_end), yrange = c(y_start, y_end))
    found_blocks[[length(found_blocks) + 1]] <- new_found_win
    
    grid_mat[r_final:(r_final + h_final - 1), c_final:(c_final + w_final - 1)] <- 0
    unassigned_ppp <- unassigned_ppp[!spatstat.geom::inside.owin(
      unassigned_ppp$x,
      unassigned_ppp$y,
      new_found_win
    )]
  }
  
  if (length(found_blocks) > 0) {
    if (verbose) cat(sprintf("\nIntegrating %d found shapes...\n", length(found_blocks)))
    
    for (found_win in found_blocks) {
      contained_ids <- unassigned_ids_current_phase[sapply(grid[unassigned_ids_current_phase], function(b) {
        small_win <- spatstat.geom::Window(b)
        all(small_win$xrange >= found_win$xrange[1] - base_tol &
              small_win$xrange <= found_win$xrange[2] + base_tol) &&
          all(small_win$yrange >= found_win$yrange[1] - base_tol &
                small_win$yrange <= found_win$yrange[2] + base_tol)
      })]
      
      if (length(contained_ids) > 0) {
        sub_blocks <- grid[contained_ids]
        x_range <- range(unlist(lapply(sub_blocks, function(x) x$window$xrange)))
        y_range <- range(unlist(lapply(sub_blocks, function(x) x$window$yrange)))
        combined_window <- spatstat.geom::owin(xrange = x_range, yrange = y_range)
        
        new_ppp <- do.call(
          spatstat.geom::superimpose,
          c(sub_blocks, list(W = combined_window))
        )
        
        success_blocks[[length(success_blocks) + 1]] <- list(ppp = new_ppp, ids = contained_ids)
        assigned[contained_ids] <- TRUE
      }
    }
  }
  
  return(list(success_blocks = success_blocks, assigned = assigned))
}





# pp_list: A list of selected block-level `ppp` objects returned by
#   `adaptive_spatial_blocking()`. Their marks are assumed to be standardized as
#   "target" and "background" for univariate analysis, or "target1", "target2",
#   and "background" for bivariate analysis.
# rratio: ratio value(s) for K-function evaluation.
#   The evaluation radius r is defined as rratio times the shorter side of each block.
# bivariate: Logical; whether to run the bivariate B-KAMP test.
bkamp_test <- function(pp_list, rratio = 0.2, bivariate = FALSE) {
  
  type1_mark = "target1"
  type2_mark = "target2"
  target_mark = "target"
  
  result_block <- lapply(seq_along(pp_list), function(index) {
    pp <- pp_list[[index]]
    win <- spatstat.geom::Window(pp)
    min_len <- min(diff(win$xrange), diff(win$yrange))
    rvec <- rratio * min_len
    n <- spatstat.geom::npoints(pp)
    
    mk <- as.character(spatstat.geom::marks(pp))
    
    if (bivariate) {
      abundance <- mean(mk %in% c(type1_mark, type2_mark))
    } else {
      abundance <- mean(mk %in% c(target_mark))
    }
    
    kamp_test(
      ppp_obj = pp,
      rvec = rvec,
      bivariate = bivariate,
      target_mark = target_mark,
      type1_mark = type1_mark,
      type2_mark = type2_mark
    ) %>%
      dplyr::mutate(
        ratio = rratio,
        block_id = index,
        n = n,
        abundance = abundance
      )
  })
  
  result_block <- dplyr::bind_rows(result_block)
  
  result_block %>%
    dplyr::group_by(ratio) %>%
    dplyr::summarise(
      combined = list(combine_stouffer(Z, weights = n / abundance)),
      b = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      statistic = purrr::map_dbl(combined, "statistic"),
      p_value = purrr::map_dbl(combined, "p_value")
    ) %>%
    dplyr::select(ratio, statistic, p_value, b)
}

# KAMP framework (Wrobel and Song, 2024)
# rvec: Numeric vector of absolute radii.
kamp_test <- function(ppp_obj, rvec, bivariate = FALSE,
                      target_mark = "target",
                      type1_mark = "target1",
                      type2_mark = "target2") {
  
  mk <- as.character(spatstat.geom::marks(ppp_obj))
  
  if (bivariate) {
    mk[!mk %in% c(type1_mark, type2_mark)] <- "background"
    mk[mk == type1_mark] <- "target1"
    mk[mk == type2_mark] <- "target2"
  } else {
    mk[mk != target_mark] <- "background"
    mk[mk == target_mark] <- "target"
  }
  
  spatstat.geom::marks(ppp_obj) <- mk
  npts  <- spatstat.geom::npoints(ppp_obj)
  areaW <- spatstat.geom:::area(spatstat.geom::Window(ppp_obj))
  
  pp_df <- as.data.frame(ppp_obj)
  distmat <- as.matrix(stats::dist(as.matrix(dplyr::select(pp_df, x, y))))
  diag(distmat) <- Inf
  
  edge_mat <- spatstat.explore::edge.Trans(ppp_obj)
  
  if (bivariate) {
    idx1 <- which(ppp_obj$marks == "target1")
    idx2 <- which(ppp_obj$marks == "target2")
    m1 <- length(idx1)
    m2 <- length(idx2)
    
    f1 <- m1 * m2 / (npts * (npts - 1))
    f2 <- f1 * (m1 + m2 - 2) / (npts - 2)
    f3 <- f1 * (m1 - 1) * (m2 - 1) / ((npts - 2) * (npts - 3))
  } else {
    idx <- which(ppp_obj$marks == "target")
    m <- length(idx)
    
    f1 <- m * (m - 1) / (npts * (npts - 1))
    f2 <- f1 * (m - 2) / (npts - 2)
    f3 <- f2 * (m - 3) / (npts - 3)
  }
  
  purrr::map_dfr(rvec, function(r) {
    Wr <- (distmat <= r) * edge_mat
    R0 <- sum(Wr)
    R1 <- sum(Wr^2)
    R2 <- sum(rowSums(Wr)^2) - R1
    R3 <- R0^2 - 2 * R1 - 4 * R2
    
    if (bivariate) {
      Kmat <- Wr[idx1, idx2, drop = FALSE]
      K <- areaW * sum(Kmat) / (m1 * m2)
      mu_K <- areaW * R0 / (npts * (npts - 1))
      var_K <- areaW^2 * (R1 * f1 + R2 * f2 + R3 * f3) /
        (m1^2 * m2^2) - mu_K^2
    } else {
      Kmat <- Wr[idx, idx, drop = FALSE]
      K <- areaW * sum(Kmat) / (m * (m - 1))
      mu_K <- areaW * R0 / (npts * (npts - 1))
      var_K <- areaW^2 * (2 * R1 * f1 + 4 * R2 * f2 + R3 * f3) /
        (m^2 * (m - 1)^2) - mu_K^2
    }
    
    Z_k <- (K - mu_K) / sqrt(var_K)
    pval <- stats::pnorm(-Z_k)
    
    tibble::tibble(
      r = r,
      khat = K,
      expectation = mu_K,
      var = var_K,
      Z = Z_k,
      pvalue = min(1, pval)
    )
  })
}


combine_stouffer <- function(z_values, weights) {
  if (length(z_values) == 0) {
    return(list(statistic = NA_real_, p_value = NA_real_))
  }
  
  statistic <- sum(weights * z_values) / sqrt(sum(weights^2))
  
  list(
    statistic = statistic,
    p_value = stats::pnorm(statistic, lower.tail = FALSE)
  )
}
