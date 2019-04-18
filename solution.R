################################################################################
# solution.R
# 
# Jeff Stevens
#
#
# Description
# -----------
#
# Functions for working with ompr solution objects.
################################################################################


source("sets.R")


# Variable description functions ------------------------------------------
#
# ompr reports solutions according to their numeric indices.  These are hard to
# keep track of, so here I translate back to the symbolic names.
# 
# 'sets' and 'vars' are convenience functions for describing the symbolic
# structure of the model variables.


# Convenience function for creating lists of sets for 'vars'
sets <- function(...) {
  # ... :   [Set name] = [Set]
  #
  # The LHS gives the name as it is to appear in the results
  
  args <- list(...)
  arg_names <- names(args)
  # Note that list(...) coerces symbolic arguments into strings...
  
  # Handle the NULL case up-front, for when the variable is unindexed.
  if ( identical(length(args), 0L) ) {
    # No arguments were specified
    return(NULL)
  }
  
  # Check that all elements of 'args' are named
  stopifnot( identical(length(arg_names), length(args)) )
  
  # Allow only integer vectors and sets as arguments
  lapply(args, function(arg) {
    stopifnot( identical(class(arg), "integer") || "set" %in% class(arg) )
  })
  
  # This list should be good as-is...Create a class out of this for future use:
  sets_bundle <- structure(args, class = "sets_bundle")
  
  return(sets_bundle)
}


# Describe the variables of a model
vars <- function(...) {
  # ... :   [variable name] = list()
  #
  # The LHS gives the variable name as used in the model
  
  args <- list(...)
  arg_names <- names(args)
  
  
  # Handle the NULL case up-front.
  # This is just to handle the pathological case of there being no variables in
  # the model, or when all variables are non-indexed.
  if ( is.null(args) ) {
    # No arguments were specified
    return(NULL)
  }
  
  
  # Check that all elements of 'args' are named
  stopifnot( identical(length(arg_names), length(args)) )
  
  # Check the argument values
  lapply(args, function(arg) {
     stopifnot( is.null(arg) || "sets_bundle" %in% class(arg) )
  })
  
  vars_bundle <- structure(args, class = "vars_bundle")
  
  return(vars_bundle)
}




# Solution extraction -----------------------------------------------------


# This extracts variable values from a solution, and translates the indices to
# the original symbolic form.
#
# 'vars_list' is the output of 'vars'; note that only solutions for the
# specified variables will be returned.
extract_solution <- function(model, solution, vars_list = NULL) {
  
  # If vars_list is NULL, then look for any non-indexed variables
  if ( is.null(vars_list) ) {
    stop( "Not implemented yet")
  }
  
  var_names <- names(vars_list)
  
  # Get the solution for a variable, and apply symbolic name indices
  get_sol <- function(vn) {
    # vn: The variable name
    
    # Get the sets bundle
    sets_list <- vars_list[[vn]]
    set_names <- names(sets_list)  # This can also be NULL
    
    # Create the ompr indexed variable descriptor expression
    var_expr <- if ( is.null(sets_list) ) {
      vn
    } else {
      paste0(vn, "[", paste(set_names, collapse = ","), "]")
    }
    
    # Get the solution
    sol <- get_solution_(solution, var_expr, type = "primal")
    ## (Should wrap this in a tryCatch statement...)
    
    # Drop the "variable" column
    if ( !is.null(sets_list) ) {
      sol[["variable"]] <- NULL
    }
   
    # Finally, translate the solution indices to their symbolic values
    #
    # The following replacement strategy isn't great practice, but is simple to
    # understand, and isn't too slow for the small number of columns used here
    if ( !is.null(set_names) ) {
      for (sn in set_names) {
        set <- sets_list[[sn]]
        
        # Only replace for symbolic sets; no need to replace integer indices
        if ( identical(class(set), "set") ) {
          # Cast this as an ordered factor for easy sorting
          keys <- names(elements(sol[[sn]], set))
          sol[[sn]] <- factor(keys, levels = names(sort(set)))

        }
      }
    }
    
    return(sol)
  }
  
  sols <- lapply( var_names, get_sol)
  sols <- setNames(sols, var_names)
  
  return(sols)
}
