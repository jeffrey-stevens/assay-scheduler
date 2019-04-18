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


# Get a subset of a set, based on a subset of keys
new_subset <- function(keys, set) {
  
  stopifnot( identical(class(set), "set") )
  
  subset <- structure( set[keys], class = "set", superset = set )
  
  return(subset)
}



# Element vectors ---------------------------------------------------------


# Create a vector of elements (order matters; replicates are permitted)
elements <- function(vals, set) {
  # vals:  A list of keys, or a list of indices
  
  if ( is.numeric(vals)) {
    # Check that all values are in the set
    stopifnot( all(vals %in% set) )
    
    # Match these values to the elements of the set
    idx <- match(vals, set)
    vec <- set[idx]
    
  } else if (is.character(vals)) {
    # Check that all the values are elements of the set
    stopifnot( all(vals %in% names(set)) )
    
    vec <- set[vals]
  
  } else {
    stop("'vals' must be an integer or character vector.")
  }
  
  obj <- structure(vec, class = "elements", set = set)
  
  return(obj)
}


# Unfortunately [] drops attributes...
`[.elements` <- function(x, i) {
  
  obj <- structure(NextMethod(), class = "elements", set = attr(x, "set"))
  
  return(obj)
}



# Map elements vectors from one set to another
#
# This defines the set mapping; it returns a function for mapping an elements
# vector on one set to the associated elements vector on the other set.
map_sets <- function(set1, set2, pairs) {
  # pairs: Pairs of keys, as a named character vector
  
  # Check the inputs
  
  # set1 and set2 must both be sets
  stopifnot( identical(class(set1), "set") && identical(class(set2), "set") )
  
  # 'pairs' must be a character vector of the same length as set1
  stopifnot( is.character(pairs) && identical(length(pairs), length(set1)) )
  
  # The names of 'pairs' must match the keys of set1, with no missing names and
  # no duplicates
  stopifnot( identical(sort(names(pairs)), sort(names(set1))) )
  
  # The values of 'pairs' are a subset of the keys of set2
  stopifnot( all(pairs %in% names(set2)) )
  
  # Todo:  Include error messages
  
  # Return a function for the association
  map <- function(vals) {
    
    # Coerce vals to an elements object, for simplicity
    keys_in <- names(elements(vals, set1))
    
    # Now map the keys of set1 to the keys of set2
    keys_out <- pairs[keys_in]
    
    # Save this as an elements vector, for simplicity
    keys_out <- elements(keys_out, set2)
    
    return(keys_out)
  }
  
  return(map)
}




# Parameters --------------------------------------------------------------
#
# Parameters in AMPL are basically lists of name-value pairs over AMPL sets.
# This is a simple implementation of this concept.
#
# Parameters can't be defined by name-value pairs like R vectors, since they
# ompr indexes by value, rather than by name...This makes the implementation
# more complex than with the "sets" defined above:

new_param <- function(keyvals, set) {
  
  # "keyvals" is a numeric vector of name-value pairs.
  
  stopifnot( is.numeric(keyvals) )
  
  # Test that all keys are in the set
  #
  # This tests that all keys are present, and that there are no duplicates:
  idx <- match(names(keyvals), names(set)) 
  stopifnot( (identical(sort(idx), seq_along(set)) ) )
  
  # Rearrange this in "set" order
  v <- keyvals[names(set)]
  
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

		# Get the parameter values
    value <- NextMethod()
    
	} else if (is.integer(i)) {
		# This is a vector of set variable indices
		set <- attr(x, "set")
		
		# Find the keys corresponding to these set values
		idx <- match(i, set)
		i <- names(set[idx])
		
		value <- NextMethod()

	} else {
		# Invalid input
		stop("Invalid indexing type.")
	}

	return(value)
}


# `[<-.param` <- function(x, i, value) {
#   
#   stopifnot( identical(length(i), length(value)))
#   
#   if (is.character(i)) {
# 		# Check that it's really a key
# 		stopifnot( all(i %in% names(x)) )
#     
#     # Set the parameter value
#     value <- NextMethod()
#     
#   } else if (is.integer(i)) {
# 		# This is a vector of set variable indices
# 		set <- attr(x, "set")
# 		
# 		# Need to keep this as i for NextMethod...
# 		i <- get_keys(i, set)
# 
# 		value <- NextMethod()
# 
# 	} else {
# 		# Invalid input
# 		stop("Invalid indexing type.")
# 	}
# 
# 	return(value)
# }
