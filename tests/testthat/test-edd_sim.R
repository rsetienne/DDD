testthat::test_that("per species rates check works", {
  testthat::expect_error(
    edd_pars_check(
      pars = c(-0.5, 0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    )
  )
  
  testthat::expect_error(
    edd_pars_check(
      pars = c(-0.5, -0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    )
  )
  
  testthat::expect_error(
    edd_pars_check(
      pars = c(0.5, -0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    )
  )
})

testthat::test_that("model and parameters match check works", {
  testthat::expect_error(
    edd_pars_check(
      pars = c(0.5, 0.1, -0.001, -0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    ),
    "this model requires six parameters"
  )
  
  testthat::expect_error(
    edd_pars_check(
      pars = c(0.5, 0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsce2",
      metric = "ed",
      offset = "none"
    ),
    "this model requires four parameters"
  )
})

testthat::test_that("metric and offset match check works", {
  testthat::expect_error(
    edd_pars_check(
      pars = c(0.5, 0.1, -0.001, -0.001),
      age = 3,
      model = "dsce2",
      metric = "ed",
      offset = "simtime"
    ),
    "only pd metric has offset methods"
  )
  
  testthat::expect_error(
    edd_pars_check(
      pars = c(0.5, 0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "spcount"
    ),
    "only pd metric has offset methods"
  )
})

testthat::test_that("Can output simulation message", {
  testthat::expect_message(
    edd_message_info(
      pars = c(0.5, 0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "spcount"
    )
  )
})

testthat::test_that("Can update rates", {
  lamu_list <- edd_update_lamu(
    ed = c(1,1,1,1),
    ed_max = c(5,5,5,5),
    params = c(4, 0.5, 0.1, -0.01, -0.01),
    model = "dsce2"
  )
  
  testthat::expect_type(
    lamu_list,
    "list"
  )
  
  testthat::expect_equal(
    length(lamu_list$newlas),
    length(lamu_list$newmus)
  )
  
  testthat::expect_equal(
    length(lamu_list$newlas),
    4
  )
  
  testthat::expect_equal(
    length(lamu_list$newmus),
    4
  )
  
  testthat::expect_false(
    NA %in% lamu_list$newmus
  )
  
  testthat::expect_false(
    NA %in% lamu_list$newlas
  )
  
  testthat::expect_false(
    NaN %in% lamu_list$newmus
  )
  
  testthat::expect_false(
    NaN %in% lamu_list$newmus
  )
  
  lamu_list <- edd_update_lamu(
    ed = c(1,1,1,1),
    ed_max = c(5,5,5,5),
    params = c(4, 0.5, 0.1, -0.01, -0.01, 0.01, 0.01),
    model = "dsde2"
  )
  
  testthat::expect_type(
    lamu_list,
    "list"
  )
  
  testthat::expect_equal(
    length(lamu_list$newlas),
    length(lamu_list$newmus)
  )
  
  testthat::expect_equal(
    length(lamu_list$newlas),
    4
  )
  
  testthat::expect_equal(
    length(lamu_list$newmus),
    4
  )
  
  testthat::expect_equal(
    NA %in% lamu_list$newmus,
    FALSE
  )
  
  testthat::expect_equal(
    NA %in% lamu_list$newlas,
    FALSE
  )
  
  testthat::expect_equal(
    NaN %in% lamu_list$newmus,
    FALSE
  )
  
  testthat::expect_equal(
    NaN %in% lamu_list$newmus,
    FALSE
  )
})

testthat::test_that("Can calculate maximum ED", {
  l_table <- structure(c(10, 10, 8.65720237313422, 7.52112256662757, 4.94259451955203, 
                        4.02141458501139, 3.18867978010257, 3.0414300418969, 2.53271705998738, 
                        1.84976965035132, 1.31034753274349, 1.26110486362692, 0.806069840362717, 
                        0.592701463871053, 0.53552669331763, 0, -1, 2, 3, 3, -1, 3, 4, 
                        -6, 5, 5, 10, 8, 12, 13, -1, 2, 3, 4, 5, -6, 7, 8, -9, 10, 11, 
                        12, 13, 14, 15, -1, -1, -1, -1, 1.26868367092781, -1, -1, -1, 
                        -1, 0.0527305518860519, -1, -1, -1, -1, -1), .Dim = c(15L, 4L
                        ))
  t <- 10
  l_table[, 1] <- t - l_table[, 1]
  num <- nrow(l_table[l_table[, 4] == -1, ])
  metric <- "pd"
  offset <- "none"
  converter <- "cpp"
  
  edmax <- edd_get_edmax(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    edmax,
    "double"
  )
  
  testthat::expect_length(
    edmax,
    num
  )
  
  offset <- "simtime"
  edmax_simtime <- edd_get_edmax(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    edmax_simtime,
    "double"
  )
  
  testthat::expect_length(
    edmax_simtime,
    num
  )
  
  offset <- "spcount"
  edmax_spcount <- edd_get_edmax(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    edmax_spcount,
    "double"
  )
  
  testthat::expect_length(
    edmax_spcount,
    num
  )
  
  offset <- "both"
  edmax_both <- edd_get_edmax(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    edmax_both,
    "double"
  )
  
  testthat::expect_length(
    edmax_both,
    num
  )
})

testthat::test_that("Can calculate ED", {
  l_table <- structure(c(10, 10, 8.65720237313422, 7.52112256662757, 4.94259451955203, 
                         4.02141458501139, 3.18867978010257, 3.0414300418969, 2.53271705998738, 
                         1.84976965035132, 1.31034753274349, 1.26110486362692, 0.806069840362717, 
                         0.592701463871053, 0.53552669331763, 0, -1, 2, 3, 3, -1, 3, 4, 
                         -6, 5, 5, 10, 8, 12, 13, -1, 2, 3, 4, 5, -6, 7, 8, -9, 10, 11, 
                         12, 13, 14, 15, -1, -1, -1, -1, 1.26868367092781, -1, -1, -1, 
                         -1, 0.0527305518860519, -1, -1, -1, -1, -1), .Dim = c(15L, 4L
                         ))
  t <- 10
  l_table[, 1] <- t - l_table[, 1]
  num <- nrow(l_table[l_table[, 4] == -1, ])
  metric <- "pd"
  offset <- "none"
  converter <- "cpp"
  
  ed <- edd_get_ed(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    ed,
    "double"
  )
  
  testthat::expect_length(
    ed,
    num
  )
  
  offset <- "simtime"
  ed_simtime <- edd_get_ed(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    ed_simtime,
    "double"
  )
  
  testthat::expect_length(
    ed_simtime,
    num
  )
  
  offset <- "spcount"
  ed_spcount <- edd_get_ed(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    ed_spcount,
    "double"
  )
  
  testthat::expect_length(
    ed_spcount,
    num
  )
  
  offset <- "both"
  ed_both <- edd_get_ed(
    num = num,
    l_table = l_table,
    t = t,
    metric = metric,
    offset = offset,
    converter <- converter
  )
  
  testthat::expect_type(
    ed_both,
    "double"
  )
  
  testthat::expect_length(
    ed_both,
    num
  )
})

testthat::test_that("Can sum rates", {
  rates <- edd_sum_rates(rep(0.5, 20), rep(0.1, 20))
  
  testthat::expect_type(
    rates,
    "double"
  )
  
  testthat::expect_length(
    rates,
    1
  )
})