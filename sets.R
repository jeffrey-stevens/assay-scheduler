################################################################################
# sets.R
# 
# Jeff Stevens
#
#
# Description
# -----------
#
# This is an attempt to implement AMPL-style indexing sets and indexed
# parameters.
#
# ompr indexes variables and parameters by integer indices.  Keeping track of
# what these indices mean becomes difficult quickly, and makes it hard to
# interpret a solution.  This provides a way of giving symbolic indices, for
# both model creation and for interpretation of the solution.
################################################################################



# Indexing sets -----------------------------------------------------------
#
# Indexing sets in AMPL provide symbolic enumerations of variable and parameter
# indices.  This is an attempt at providing some symbolic indexing with simple
# subsetting characteristics.


new_set <- function(keys) {
  # keys:  A character vector, or a name-value character vector
  
  # Check that keys is a character or integer vector
  stopifnot( is.character(keys) || is.integer(keys) )
  
  if (is.character(keys)) {
    # 'keys' is a character vector; generate a sequence of index values
    
    # There shouldn't be any duplicate keys
    if ( !identical(anyDuplicated(keys), 0L) ) {
      stop("'key' musn't contain duplicate values.")
    }
    
    values <- seq_along(keys)
    set <- structure( setNames(seq_along(keys), keys),
                      class = "set", superset = NULL )
    # Setting an attribute to NULL has no effect (the attribute isn't stored),
    # but it maintains consistency with the subset structure.  Unfortunately
    # there doesn't appear to be a way of setting a NULL attribute...
    
  } else if (is.integer(keys)) {
    
    # Check that all keys are named
    key_names <- names(keys)
    if ( !identical(length(key_names), length(keys)) ) {
      stop("Some names are missing from 'keys'.")
    }
    
    # Check that the names and values are unique
    if ( !identical(anyDuplicated(key_names), 0L) ) {
      stop("'key' musn't contain duplicate names.")
    }
    if ( !identical(anyDuplicated(keys), 0L) ) {
      stop("'key' musn't contain duplicate values.")
    }
    
    set <- structure(keys, class = "set", superset = NULL)
  } else {
    
    stop("'keys' must either be a character vector, or a named integer vector.")
  }
  
  
  return(set)
}


# Get the keys corresponding to a vector of variable indices (columns)
get_keys <- function(set, vals = set) {
  
  # Check that the indices in 'vals' are within 'set'
  stopifnot( all(vals %in% set) )
  
  # Get the indices of the set items in vals 
  idx <- match(vals, set)
  # Note that there shouldn't be any NAs in here...
  
  # The keys are the names of the corresponding subset
  keys <- names(set[idx])
  
  return(keys)
}


# Get a subset of a set, based on a subset of keys
new_subset <- function(set, keys) {
  
  stopifnot( identical(class(set), "set") )
  
  subset <- structure( set[keys], class = "set", superset = set )
  
  return(subset)
}


# Indexing should return a subset
`[.set` <- function(x, i) {
  
  subset <- structure(NextMethod(), class = "set", superset = x)
  
  return(subset)
}


# Create a vector of elements (order matters; replicates are permitted)
element_vector <- function(vals, set) {
  # vals:  A list of keys, or a list of indices
  
  if ( is.numeric(vals)) {
    # Get the keys corresponding to these values
    keys <- get_keys(set, vals)
    
  } else if (is.character(vals)) {
    # Check that all the values are elements of the set
    stopifnot( all(vals %in% get_keys(set)) )
    
    keys <- vals
  
  } else {
    stop("'vals' must be an integer or character vector.")
  }
  
  # Get the elements, clearing the attributes (except for "names")
  vec <- unclass(set[keys])
  attr(vec, "superset") <- NULL
  
  obj <- structure(vec, class = "element_vec", set = set)
  
  return(obj)
}


# Unfortunately [] drops attributes...
`[.element_vec` <- function(x, i) {
  
  obj <- structure(NextMethod(), class = "element_vec", set = attr(x, "set"))
  
  return(obj)
}


# Get the next element of a set
shift_set <- function(set, within = TRUE) {
  # within:  Keep within the subset, or shift into the superset?
  
  if (within) {
    # Create a vector of indices
    idx <- seq_along(set)
    
    # Shift them forward by one, truncating the last element
    shifted_idx <- (idx + 1)[-length(idx)]
    
    # Return the new values as a subset
    shifted_set <- new_subset(set, get_keys(set)[shifted_idx])
    
  } else {
    # Shifting out of the set is more involved...
    superset <- attr(set, "superset")
    N <- length(superset)
    
    idx <- match(set, superset)
    shifted_idx <- idx + 1
    n <- length(shifted_idx)
    
    # Drop the last index if it's out of range
    if (shifted_idx[[n]] > N) {
      shifted_idx <- shifted_idx[-n]
    }
    
    # Return the new values as a subset of the superset
    shifted_set <- new_subset(superset, get_keys(superset)[shifted_idx])
  }
  
  return(shifted_set)
}


# Shift a vector (ompr strips attributes, unfortunately:
shift_vec <- function(vec, set) {
  
  keys <- get_keys(set, vec)
  subset <- new_subset(set, keys)
  
  shifted_subset <- shift_set(subset, within = FALSE)
  
  # Strip this of any structure for ompr
  return( as.vector(shifted_subset) )
  
}


# Map between two sets
map_sets <- function(set1, set2, pairs) {
  # pairs: Pairs of keys, as a named character vector
  
  # Check the inputs
  
  # set1 and set2 must both be sets
  stopifnot( identical(class(set1), "set") && identical(class(set2), "set") )
  
  # 'pairs' must be a character vector of the same length as set1
  stopifnot( is.character(pairs) && identical(length(pairs), length(set1)) )
  
  # The names of 'pairs' must match the keys of set1, with no missing names and
  # no duplicates
  stopifnot( identical(sort(names(pairs)), sort(get_keys(set1))) )
  
  # The values of 'pairs' are a subset of the keys of set2
  stopifnot( all(pairs %in% get_keys(set2)) )
  
  # Todo:  Include error messages
  
  # Return a function for the association
  map <- function(vals1) {
    
    # Check that vals is a subset of set1
    stopifnot( all(vals1 %in% set1) )
    
    # First, get the keys of set1 associated with vals
    keys1 <- get_keys(set1, vals1)
    
    # Now map the keys of set1 to the keys of set2
    keys2 <- pairs[keys1]
    
    # Return the result as a vector
    vals2 <- as.vector(set2[keys2])
    
    return(vals2)
  }
  
  return(map)
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
