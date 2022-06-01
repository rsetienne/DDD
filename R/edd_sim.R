#' @name create_edd_pars
#' @title Create a list of parameter sets for edd simulation
#' @description Function to conveniently create a list containing all possible
#' parameter combinations for edd simulation
#' @param save_file Logical, decides whether to save the created list to file
#' @param ... parameters that form a parameter space
#' @return Returns a list containing the parameter sets
#' @author Tianjian Qin
#' @importFrom eve edd_combo_maker
#' @export create_edd_pars
create_edd_pars <- function(save_file = FALSE, ...) {
  return(edd_combo_maker(save_file, ...))
}



#' @name edd_simulation
#' @title Start an edd simulation
#' @description Function to start an edd simulation, of which the evolution
#' process is dependent on phylogenetic relationship among the community
#' @param combo A list containing parameter sets to be used to start a simulation
#' @param history Logical, decides whether to record historical states for all
#' of the lineages, including the transitions of the rates and the evolutionary
#' relationship metrics
#' @param verbose Logical, decides whether to print simulation details of each
#' time step during simulation, for debugging purposes
#' @param nrep Number of replications
#' @return Returns a large list containing all the results
#' @author Tianjian Qin
#' @export edd_simulation
#' @importFrom eve edd_sim_rep
edd_simulation <-
  function(combo = NULL,
           history = FALSE,
           verbose = FALSE,
           nrep = 5) {
    lapply(
      combo,
      edd_sim_rep,
      history = history,
      verbose = verbose,
      nrep = nrep
    )
  }