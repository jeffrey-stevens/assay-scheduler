################################################################################
# scheduler.R
# 
# Jeff Stevens
#
#
# Description
# -----------
#
# This is a recreation of a MIP scheduler that I had written in opmr/GLPK to
# help determine optimal incubation times and overall throughput for automating
# a magnetic bead assay.
# 
# Note that this is an exponential solution to the problem, so can be terribly
# slow for more than a few plates.  Still, it served us well in finding optimal
# throughput schedules for automating our assay on a liquid handling system.
################################################################################


source("sets.R")


library(dplyr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)


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


WASH_STEPS <- new_subset(STEPS, c("WashBeads", "WashSample",
                                  "WashPrimary", "WashSecondary"))

INCUBATION_STEPS <- new_subset(STEPS, c("IncubateSample", "IncubatePrimary",
                                      "IncubateSecondary"))

PROCESSING_STEPS <- new_subset(STEPS, 
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
    ##
    ## Keep the incubation times consistent between plates
    add_variable(inc_dur[step],
                 step = INCUBATE_STEPS,
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

  ## Duration constraint
  
  durations <- new_param(PROCESSING_STEPS,
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
  
  model <-
    model %>%
    add_constraint( begin_time[plate, step] + durations[step] == end_time[plate, step],
                    plate = PLATES, step = PROCESSING_STEPS )

  ## Step ordering constraints
  
  model <-
    model %>%
    ## Order the processing steps
    add_constraint( end_time[plate, step1] <= begin_time[plate, step2],
                    plate = PLATES,
                    step1 = PROCESSING_STEPS[-length(PROCESSING_STEPS)],
                    step2 = PROCESSING_STEPS[-1],
                    # ompr removes all attributes, so must use a vectorized version
                    step2 == shift_vec(step1, PROCESSING_STEPS) )
  
  ## Incubation definition constraints
  
  ### Incubations need an equality constraint, since they're defined as the gap between
  ### two processing steps...
  
  ### Creating set maps makes things easy to comprehend...
  
  ### Incubations start with the addition of reagent, and end with washing the plate
  inc_start_map <- map_sets(INCUBATION_STEPS, PROCESSING_STEPS,
                            c("IncubateSample" = "AddSample",
                              "IncubatePrimary" = "AddPrimary",
                              "IncubateSecondary" = "AddSecondary"))
  inc_end_map <- map_sets(INCUBATION_STEPS, PROCESSING_STEPS,
                          c("IncubateSample" = "WashSample",
                            "IncubatePrimary" = "WashPrimary",
                            "IncubateSecondary" = "WashSecondary"))
  
  model <-
    model %>%
    add_constraint( inc_dur[step] == begin_time[plate, step2] - end_time[plate, step1],
                    plate = PLATES,
                    step = INCUBATION_STEPS, step1 = PROCESSING_STEPS, step2 = PROCESSING_STEPS,
                    step1 == inc_start_map(step), step2 == inc_end_map(step))
  
  ### For now, require that the incubations are the same for all plates
  inc_dur_min <- new_param(INCUBATION_STEPS, c("IncubateSample" = params$sample_inc_min,
                                               "IncubatePrimary" = params$primary_inc_min,
                                               "IncubateSecondary" = params$secondary_inc_min))
  inc_dur_max <- new_param(INCUBATION_STEPS, c("IncubateSample" = params$sample_inc_max,
                                               "IncubatePrimary" = params$primary_inc_max,
                                               "IncubateSecondary" = params$secondary_inc_max))

  
  ## Assay and run constraints



  # ----- Objective -----

  
  return(model)
}



solve_schedule <- function() {
}

