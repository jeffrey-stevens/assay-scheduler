library(dplyr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)



# Indexing sets and indexed parameters -----------------------------------------------------------


# ----- Implementation -----

# OMPR doesn't implement the concept of indexing sets.
# This implements something close to it

new_set <- function(keys) {
  # keys:  A character vector
  
  set <- setNames(seq_along(keys), keys)
  
  return(set)
}


# Get the keys corresponding to a vector of variable indices (columns)
get_keys <- function(set, varidx = set) {
  
  # Get the indices of the set items in varidx 
  idx <- which(set %in% varidx)
  
  # The keys are the names of the corresponding subset
  keys <- names(set[idx])
  
  return(keys)
}

# Get a subset of a set, based on a subset of keys
get_subset <- function(set, keys) {
  set[keys]
}


# Parameters are more difficult to implement since they may be defined over
# subsets, yet ompr prefers to index by position...

new_param <- function(set, values) {
  
  # "values" can be missing, NULL, or a vector of key-value pairs,
  # where the keys are among the keys of the set.  Set any missing
  # values to NA.
  
  if ( missing(values) || is.null(values) ) {
    # Create an empty vector of NAs (inefficient, but clean) 
    v <- setNames( rep(NA_real_, length(set)), get_keys(set) )
  
  } else if ( identical(length(names(values)), length(set)) ) {
    
    if ( ! (all(names(values) %in% get_keys(set))) ) {
      stop("All names of 'values' must be in 'set'.")
    }
    
    # Rearrange this in "set" order
    v <- values[get_keys(set)]
    
    # Of course this gets messy if there are duplicate names...
  }
  
  else {
    stop("'values' must either be missing, NULL, or a named vector with names from 'set'.")
  }
  
	# This just creates a vector of values indexed by set element key.
	# Making this into a class and attaching the set as an attribute allows
	# `[` to be overridden...
	param <- structure(v, class = "param", set = set)

	return(param)
}



`[.param` <- function(x, i) {
	
	# If i is a character vector, then assume that's a set element key:
	if (is.character(i)) {
		# Check that it's really a key
		stopifnot( all(i %in% names(x)) )

		# Get the parameter value
	  # x needs to be unclassed to prevent infinite recursion...
		value <- NextMethod()

	} else if (is.integer(i)) {
		# This is a vector of set variable indices
		set <- attr(x, "set")
		
		# Need to keep this as i for NextMethod...
		i <- get_keys(set, i)

		value <- NextMethod()

	} else {
		# Invalid input
		stop("Invalid indexing type.")
	}

	return(value)
}


`[<-.param` <- function(x, i, value) {
  
  stopifnot( identical(length(i), length(value)))
  
  if (is.character(i)) {
		# Check that it's really a key
		stopifnot( all(i %in% names(x)) )
    
    # Set the parameter value
    value <- NextMethod()
    
  } else if (is.integer(i)) {
		# This is a vector of set variable indices
		set <- attr(x, "set")
		
		# Need to keep this as i for NextMethod...
		i <- get_keys(set, i)

		value <- NextMethod()

	} else {
		# Invalid input
		stop("Invalid indexing type.")
	}

	return(value)
}



# ----- Sets -----

STEPS <- new_set( 
  c(  "AddBeads",
      "WashBeads",
      "AddSample",
      "IncubateSample",
      "WashSample",
      "AddPrimary",
      "IncubatePrimary",
      "WashPrimary",
      "AddSecondary",
      "IncubateSecondary",
      "WashSecondary",
      "AddEnhancer",
      "Read"
  )
)


WASH_STEPS <- get_subset(STEPS, c("WashBeads", "WashSample",
                                  "WashPrimary", "WashSecondary"))

INCUBATION_STEPS <- get_subset(STEPS, c("IncubateSample", "IncubatePrimary",
                                      "IncubateSecondary"))

PROCESSING_STEPS <- get_subset(STEPS, 
                               setdiff(get_keys(STEPS), get_keys(INCUBATION_STEPS)))


# Defaults ----------------------------------------------------------------

# Timings are in *seconds*

DEFAULTS <- list(
  add_beads_dur = 120,
  add_sample_dur = 60 * 7,
  add_primary_dur = 90,
  add_secondary_dur = 90,
  add_enhancer_dur = 90,
  wash_dur = 90,
  read_dur = 60 * 40,
  sample_inc_min = 60 * 25,
  sample_inc_max = 60 * 35,
  primary_inc_min = 60 * 20,
  primary_inc_max = 60 * 35,
  secondary_inc_min = 60 * 10,
  secondary_inc_max = 60 * 35
)



# Scheduler ---------------------------------------------------------------


build_model <- function(num_plates = 1L, num_readers = 1L,
                        regular = TRUE, params = DEFAULTS) {
  
  # ----- Setup -----
  PLATES <- seq_len(num_plates)
  
  DURATIONS <- new_param(PROCESSING_STEPS,
                         c("AddBeads" = params$add_beads_dur,
                           "WashBeads" = params$wash_dur,
                           "AddSample" = params$add_sample_dur,
                           "WashSample" = params$wash_dur,
                           "AddPrimary" = params$add_primary_dur,
                           "WashPrimary" = params$wash_dur,
                           "AddSecondary" = params$add_secondary_dur,
                           "WashSecondary" = params$wash_dur,
                           "AddEnhancer" = params$add_enhancer_dur,
                           "Read" = params$read_dur) )
  
  model <- MIPModel()
  
  
  # ----- Define the variables -----
  
  model <-
    model %>%
    
    ## Step start and stop times
    add_variable(begin_time[plate, step],
                 plate = PLATES, step = STEPS,
                 type = "continuous", lb = 0) %>%
    add_variable(end_time[plate, step],
                 plate = PLATES, step = STEPS,
                 type = "continuous", lb = 0) %>%
    
    ## Incubation durations
    add_variable(dur[plate, step],
                 plate = PLATES, step = INCUBATE_STEPS,
                 type = "continuous", lb = 0) %>%
    
    ## Assay timings and durations
    add_variable(assay_begin[plate],
                 plate = PLATES,
                 type = "continuous", lb = 0) %>%
    add_variable(assay_end[plate],
                 plate = PLATES,
                 type = "continuous", lb = 0) %>%
    add_variable(assay_dur[plate],
                 plate = PLATES,
                 type = "continuous", lb = 0) %>%
    
    ## Total run end time and duration (start time is assumed to be at t = 0) %>%
    add_variable(run_end, type = "continuous", lb = 0) %>%
    add_variable(run_dur, type = "continuous", lb = 0)
  

  # ----- Constraints -----

  ## Duration constraints
  DURATIONS <- new_param(PROCESSING_STEPS)
  


  ## Step ordering constraints



  # ----- Objective -----

  
  return(model)
}



solve_schedule <- function() {
}

