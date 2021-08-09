edd_pars_check <- function(pars, age, model, metric, offset) {
  # pars range check
  if (pars[1] <= 0 | pars[2] <= 0) {
    stop("per species rates should be positive")
  }
  if (pars[2] <= 0) {
    stop("coefficient for extinction should be positive")
  }
  # pars and model match check
  if (model == "dsce2" && length(pars) != 4) {
    stop("this model requires four parameters")
  }
  if (model == "dsde2" && length(pars) != 6) {
    stop("this model requires six parameters")
  }
  # metric and offset match check
  if (metric != "pd" && offset != "none") {
    stop("only pd metric has offset methods")
  }
}

edd_update_lamu <- function(ed, ed_max, params, model) {
  num <- params[1]
  la0 <- params[2]
  mu0 <- params[3]
  beta_num <- params[4]
  beta_phi <- params[5]
  if (model == "dsce2") {
    if (length(params) != 5) {
      stop("incorrect parameter(s)")
    }
    # dependent speciation, constant extinction
    if (beta_phi < 0) {
      newlas <- pmax(0, la0 + beta_num * num + beta_phi * ed)
    } else {
      newlas <- pmax(0, la0 + beta_num * num + beta_phi * ed_max)
    }
    newmus <- rep(mu0, length(newlas))
  } else if (model == "dsde2") {
    if (length(params) != 7) {
      stop("incorrect parameter(s)")
    }
    # dependent speciation, dependent extinction
    gamma_num <- params[6]
    gamma_phi <- params[7]
    if (beta_phi < 0) {
      newlas <- pmax(0, la0 + beta_num * num + beta_phi * ed)
    } else {
      newlas <- pmax(0, la0 + beta_num * num + beta_phi * ed_max)
    }
    if (gamma_phi < 0) {
      newmus <- pmax(0, mu0 + gamma_num * num + gamma_phi * ed)
    } else {
      newmus <- pmax(0, mu0 + gamma_num * num + gamma_phi * ed_max)
    }
  }
  return(list(newlas = newlas, newmus = newmus))
}

edd_get_edmax <- function(num, l_table, t, metric, offset) {
  if (metric == "ed") {
    ed_max <- as.vector(DDD::L2ED(l_table, t))
  } else if (metric == "pd") {
    if (offset == "none") {
      ed_max <- rep(as.vector(DDD::L2Phi(l_table, t, metric)), num)
    } else if (offset == "simtime") {
      ed_max <- rep(as.vector(DDD::L2Phi(l_table, t, metric) - t), num)
    } else {
      stop("no such offset method")
    }
  } else {
    stop("no such metric")
  }
  return(ed_max)
}

edd_get_ed <- function(num, l_table, t, metric, offset) {
  if (metric == "ed") {
    ed <- as.vector(DDD::L2ED(l_table, t))
  } else if (metric == "pd") {
    if (offset == "none") {
      ed <- rep(as.vector(DDD::L2Phi(l_table, t, metric)), num)
    } else if (offset == "simtime") {
      ed <- rep(as.vector(DDD::L2Phi(l_table, t, metric) - t), num)
    } else {
      stop("no such offset method")
    }
  } else {
    stop("no such metric")
  }
  return(ed)
}

edd_sum_rates <- function(las, mus) {
  return(sum(las) + sum(mus))
}

edd_sample_event <- function(las, mus, linlist) {
  events <- 1:(2 * length(linlist))
  return(DDD::sample2(events, 1, prob = c(las, mus)))
}

#' @title Function to simulate the evolutionary distinctiveness dependent
#' diversification process
#' @description Simulating the evolutionary distinctiveness dependent
#' diversification process
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to
#' lambda (speciation rate) \cr \code{pars[2]} corresponds to mu (extinction
#' rate) \cr \code{pars[3]} corresponds to beta_num (coefficient for species
#' number effect on speciation) \cr \code{pars[4]} corresponds to beta_phi
#' (coefficient for evolutionary distinctiveness effect on speciation)
#' \cr \code{pars[5]} corresponds to gamma_num (coefficient for species number
#' effect on speciation) \cr \code{pars[6]} corresponds to gamma_phi
#' (coefficient for evolutionary distinctiveness effect on extinction)
#' @param age Sets the crown age for the simulation
#' @param model Sets the model of diversity-dependence: \cr \code{model ==
#' dsce2} : linear dependence in speciation rate with parameters beta_num and
#' beta-phi\cr \code{model == dsde2} : linear dependence
#' in both speciation rate and extinction rate with parameters beta_num,
#' beta_phi, gamma_num and gamma_phi
#' @param metric Specifies which phylogenetic diversity metric should be used
#' @param offset Specifies which method to use to offset time effect on
#' evolutionary distinctiveness
#' @return \item{ out }{ A list with the following nine elements: The first
#' element is the tree of extant species in phylo format \cr The second element
#' is the tree of all species, including extinct species, in phylo format \cr
#' The third element is a matrix of all species where \cr - the first column is
#' the time at which a species is born \cr - the second column is the label of
#' the parent of the species; positive and negative values only indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values only indicate whether the species belongs to the left or
#' right crown lineage \cr - the fourth column is the time of extinction of the
#' species. If this equals -1, then the species is still extant.\cr The fourth
#' element is the set of branching times of the tree of extant species.\cr The
#' fifth element is the lineage-through-time plot. \cr The sixth element is a
#' list of all evolutionary distinctiveness values of all lineages. \cr The
#' seventh element is a list of all the speciation rates of all lineages at all
#' the time steps. \cr The eighth element is a list of all the extinction rates
#' of all lineages at all the time steps. \cr The nineth element is a list of
#' all the lineages at all the time steps.}
#' @author Tianjian Qin, Rampal S. Etienne
#' @keywords models
#' @examples
#' edd_sim(
#'   pars = c(0.5, 0.1, -0.001, -0.001, 0.001, 0.001), age = 6,
#'   model = "dsde2", metric = "ed", offset = "none"
#' )
#' @export edd_sim
edd_sim <- function(pars,
                    age,
                    model = "dsce2",
                    metric = "ed",
                    offset = "none") {
  edd_pars_check(pars, age, model, metric, offset)

  done <- 0
  while (done == 0) {
    # initialization
    t <- rep(0, 1)
    l_table <- matrix(0, 2, 4)
    i <- 1
    t[1] <- 0
    num <- 2
    l_table[1, 1:4] <- c(0, 0, -1, -1)
    l_table[2, 1:4] <- c(0, -1, 2, -1)
    linlist <- c(-1, 2)
    new_lin <- 2
    params <- c(num, pars)
    ed <- c(0, 0)
    ed_max <- edd_get_edmax(num, l_table, age, metric, offset)
    lamu <- edd_update_lamu(ed, ed_max, params, model)

    # controlling significant digits in tibble objects
    options(pillar.sigfig = 10)

    # store EDs and lamus and associated lineages
    eds <- list(ed)
    las <- list(lamu$newlas)
    mus <- list(lamu$newmus)
    linlists <- list(linlist)

    # get time interval
    t[i + 1] <-
      t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))

    # main simulation circle
    while (t[i + 1] <= age) {
      # time step index
      i <- i + 1
      ed <- edd_get_ed(num[i - 1], l_table, t[i], metric, offset)
      lamu_real <- edd_update_lamu(ed, ed, params, model)
      event_type <- sample(c("real", "fake"),
        1,
        prob = c(
          sum(lamu_real$newlas + lamu_real$newmus),
          sum(
            lamu$newlas - lamu_real$newlas + lamu$newmus - lamu_real$newmus
          )
        )
      )
      if (event_type == "real") {
        event <-
          edd_sample_event(lamu_real$newlas, lamu_real$newmus, linlist)
        ran_lin <- c(linlist, linlist)[event]

        if (event <= length(linlist)) {
          num[i] <- num[i - 1] + 1
          new_lin <- new_lin + 1
          l_table <-
            rbind(l_table, c(t[i], ran_lin, sign(ran_lin) * new_lin, -1))
          linlist <- c(linlist, sign(ran_lin) * new_lin)
        } else {
          num[i] <- num[i - 1] - 1
          l_table[abs(ran_lin), 4] <- t[i]
          w <- which(linlist == ran_lin)
          linlist <- linlist[-w]
        }
      } else {
        num[i] <- num[i - 1]
      }

      if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
        t[i + 1] <- Inf
      } else {
        ed <- edd_get_ed(num[i], l_table, t[i], metric, offset)
        ed_max <- edd_get_edmax(num[i], l_table, age, metric, offset)
        params[1] <- num[i]
        lamu <- edd_update_lamu(ed, ed_max, params, model)

        # append to EDs and lamus
        eds <- c(eds, list(ed))
        las <- c(las, list(lamu$newlas))
        mus <- c(mus, list(lamu$newmus))
        linlists <- c(linlists, list(linlist))

        t[i + 1] <-
          t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))
      }
    }

    if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
      done <- 0
    } else {
      done <- 1
    }
  }

  tes <- DDD::L2phylo2(l_table, age, dropextinct = T)
  tas <- DDD::L2phylo2(l_table, age, dropextinct = F)
  brts <- DDD::L2brts2(l_table, age, dropextinct = T)
  nltt <-
    data.frame(
      "time" = t[-i],
      "num" = num
    )
  nltt[nrow[nltt], 1] <- age

  out <-
    list(
      tes = tes,
      tas = tas,
      l_table = l_table,
      brts = brts,
      nltt = nltt,
      eds = eds,
      las = las,
      mus = mus,
      linlists = linlists
    )

  return(out)
}