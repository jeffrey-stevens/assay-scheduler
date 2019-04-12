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




# Constants and globals ---------------------------------------------------


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



# ----- Defaults -----

# Timings are in *seconds*

DEFAULTS <- list(
  ## Time to mix the bead suspension, and transfer it into a full plate
  add_beads_dur = 120,
  
  ## Time to transfer sample from 96 test tubes into a full plate
  add_sample_dur = 60 * 7.5,
  
  ## Time to transfer the primary detector into a full plate
  add_primary_dur = 90,
  
  ## Time to transfer the secondary detector into a full plate
  add_secondary_dur = 90,
  
  ## Time to transfer enhancement buffer (for signal stabilization and
  ## enhancement) into a full plate
  add_enhancer_dur = 90,
  
  ## The duration of a wash step, assuming that the same wash protocol is used
  ## for each wash
  wash_dur = 90,
  
  ## The duration of a read step.  Our reader was excruciatingly slow, due to
  ## the unique technology involved.
  read_dur = 60 * 40,
  
  ## Bounds on the incubation time for sample incubation
  sample_inc_min = 60 * 25,
  sample_inc_max = 60 * 35,
  
  ## Bounds on the incubation time for primary detector incubation
  primary_inc_min = 60 * 20,
  primary_inc_max = 60 * 35,
  
  ## Bounds on the incubation time for secondary detector incubation
  secondary_inc_min = 60 * 10,
  secondary_inc_max = 60 * 35
  
  ## Enhancement buffer doesn't require any incubation
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
    add_variable(run_end, type = "continuous", lb = 0)
  

  # ===== Constraints =====

  
  # ----- Duration constraint ----- 
  
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

  
  # ----- Step ordering constraints -----
  
  # Explicitly map the sequence of steps
  next_step <- map_sets(PROCESSING_STEPS[-length(PROCESSING_STEPS)],
                        PROCESSING_STEPS[-1],
                        c( "AddBeads" = "WashBeads",
                           "WashBeads" = "AddSample",
                           "AddSample" = "WashSample",
                           "WashSample" = "AddPrimary",
                           "AddPrimary" = "WashPrimary",
                           "WashPrimary" = "AddSecondary",
                           "AddSecondary" = "WashSecondary",
                           "WashSecondary" = "AddEnhancer",
                           "AddEnhancer" = "Read"))
  
  model <-
    model %>%
    ## Order the processing steps
    add_constraint( end_time[plate, step1] <= begin_time[plate, step2],
                    plate = PLATES,
                    step1 = PROCESSING_STEPS[-length(PROCESSING_STEPS)],
                    step2 = PROCESSING_STEPS[-1],
                    # ompr removes all attributes, so must use a vectorized version
                    step2 == next_step(step1) )
  
  
  
  # ----- Incubation definitions -----
  
  # Incubations need an equality constraint, since they're defined as the gap between
  # two processing steps...
  
  # Creating set maps makes things easy to comprehend...
  
  # Incubations start with the addition of reagent, and end with washing the plate
  inc_start_map <- map_sets(INCUBATION_STEPS, PROCESSING_STEPS,
                            c("IncubateSample" = "AddSample",
                              "IncubatePrimary" = "AddPrimary",
                              "IncubateSecondary" = "AddSecondary"))
  inc_end_map <- map_sets(INCUBATION_STEPS, PROCESSING_STEPS,
                          c("IncubateSample" = "WashSample",
                            "IncubatePrimary" = "WashPrimary",
                            "IncubateSecondary" = "WashSecondary"))
  
  # browser()
  
  model <-
    model %>%
    add_constraint( begin_time[plate, inc_step] == end_time[plate, proc_step],
                    plate = PLATES, inc_step = INCUBATION_STEPS, proc_step = PROCESSING_STEPS,
                    proc_step == inc_start_map(inc_step) ) %>%

    add_constraint( end_time[plate, inc_step] == begin_time[plate, proc_step],
                    plate = PLATES, inc_step = INCUBATION_STEPS, proc_step = PROCESSING_STEPS,
                    proc_step == inc_end_map(inc_step) ) %>%

    add_constraint( inc_dur[step] == end_time[plate, step] - begin_time[plate, step],
                    ## All incubations are assumed to be the same (for now)
                    plate = PLATES[1], step = INCUBATION_STEPS )
  
  
  # ----- Incubation bounds constraints -----
  
  # For now, require that the incubations are the same for all plates
  inc_dur_min <- new_param(INCUBATION_STEPS, c("IncubateSample" = params$sample_inc_min,
                                               "IncubatePrimary" = params$primary_inc_min,
                                               "IncubateSecondary" = params$secondary_inc_min))
  inc_dur_max <- new_param(INCUBATION_STEPS, c("IncubateSample" = params$sample_inc_max,
                                               "IncubatePrimary" = params$primary_inc_max,
                                               "IncubateSecondary" = params$secondary_inc_max))

  model <-
    model %>%
    add_constraint( inc_dur[step] >= inc_dur_min[step],
                    step = INCUBATION_STEPS) %>%
    add_constraint( inc_dur[step] <= inc_dur_max[step],
                    step = INCUBATION_STEPS)
  
  
  # ----- Assay and run definitions -----
  
  # An assay begins with the beginning of bead suspension transfer step, and
  # ends with the end of the read (note that plate transport steps are being
  # ignored for now...)
  
  model <-
    model %>%
    add_constraint( assay_begin[plate] == begin_time[plate, step],
                    plate = PLATES, step = STEPS["AddBeads"] ) %>%
    add_constraint( assay_end[plate] == end_time[plate, step],
                    plate = PLATES, step = STEPS["Read"]) %>%
    add_constraint( assay_dur[plate] == assay_end[plate] - assay_begin[plate],
                    plate = PLATES)
  
  # A run ends with the end of the last assay
  
  model <-
    model %>%
    add_constraint( run_end == assay_end[plate],
                   plate = PLATES[length(PLATES)] )


  # ===== Objective =====
  
  # Minimize the overall run time for this number of plates
  model <- set_objective(model, run_end, sense = "min")

  
  return(model)
}



solve_schedule <- function(model) {
  
  solution <- solve_model(model, with_ROI(solver = "glpk"))
  
  return(solution)
}
