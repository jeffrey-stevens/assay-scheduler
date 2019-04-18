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
source("solution.R")


library(dplyr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)
library(dplyr)



# Constants and globals ---------------------------------------------------


# ----- Resources -----

# Define the physical resources used by the system
RESOURCES <- new_set( c( "Pipettor", "Washer", "Reader") )


# ----- Steps -----

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

WASH_STEPS <- new_subset(c("WashBeads", "WashSample",
                           "WashPrimary", "WashSecondary"),
                         STEPS)

INCUBATION_STEPS <- new_subset(c("IncubateSample", "IncubatePrimary",
                                 "IncubateSecondary"),
                               STEPS)

PROCESSING_STEPS <- new_subset( setdiff(names(STEPS), names(INCUBATION_STEPS)),
                                STEPS)


# ----- Maps -----

# Map resources 
resources_map <- map_sets( PROCESSING_STEPS, RESOURCES, 
  c(  "AddBeads"          = "Pipettor",
      "WashBeads"         = "Washer",
      "AddSample"         = "Pipettor",
      "WashSample"        = "Washer",
      "AddPrimary"        = "Pipettor",
      "WashPrimary"       = "Washer",
      "AddSecondary"      = "Pipettor",
      "WashSecondary"     = "Washer",
      "AddEnhancer"       = "Pipettor",
      "Read"              = "Reader"
  )
)


# ----- Defaults -----

# Timings are in *seconds*

# The default list of parameters for build_model
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
  secondary_inc_max = 60 * 35,
  
  ## Enhancement buffer doesn't require any incubation
  
  ## The maximum time a plate can be "exposed" after washing
  wash_exposure_max = 2 * 60,
  
  ## The maximum time beads can sit in enhancement buffer before being read
  enhancer_exposure_max = 10 * 60
  
)



# Scheduler ---------------------------------------------------------------


build_model <- function(num_plates = 1L, num_readers = 1L,
                        regular = TRUE, params = DEFAULTS) {
  # regular :  Run the plates at regular intervals?
  
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
                 step = INCUBATION_STEPS,
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
  
  durations <- new_param( c("AddBeads" = params$add_beads_dur,
                           "WashBeads" = params$wash_dur,
                           "AddSample" = params$add_sample_dur,
                           "WashSample" = params$wash_dur,
                           "AddPrimary" = params$add_primary_dur,
                           "WashPrimary" = params$wash_dur,
                           "AddSecondary" = params$add_secondary_dur,
                           "WashSecondary" = params$wash_dur,
                           "AddEnhancer" = params$add_enhancer_dur,
                           "Read" = params$read_dur),
                          PROCESSING_STEPS )
  
  model <-
    model %>%
    add_constraint( begin_time[plate, step] + durations[step] == end_time[plate, step],
                    plate = PLATES, step = PROCESSING_STEPS )

  
  # ----- Step ordering constraints -----
  
  # Explicitly map the sequence of steps
  
  pre_steps <- new_subset( c( "AddBeads",
                              "WashBeads",
                              "AddSample",
                              "WashSample",
                              "AddPrimary",
                              "WashPrimary",
                              "AddSecondary",
                              "WashSecondary",
                              "AddEnhancer"),
                           PROCESSING_STEPS )
  
  post_steps <- new_subset( c(  "WashBeads",
                                "AddSample",
                                "WashSample",
                                "AddPrimary",
                                "WashPrimary",
                                "AddSecondary",
                                "WashSecondary",
                                "AddEnhancer",
                                "Read"),
                            PROCESSING_STEPS )
  
  next_step <- map_sets(pre_steps,
                        post_steps,
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
                    step1 = pre_steps,
                    step2 = post_steps,
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
                    plate = PLATES, step = INCUBATION_STEPS )
  
  
  # ----- Incubation bounds constraints -----
  
  # For now, require that the incubations are the same for all plates
  inc_dur_min <- new_param( c("IncubateSample" = params$sample_inc_min,
                              "IncubatePrimary" = params$primary_inc_min,
                              "IncubateSecondary" = params$secondary_inc_min),
                            INCUBATION_STEPS )
  inc_dur_max <- new_param( c("IncubateSample" = params$sample_inc_max,
                              "IncubatePrimary" = params$primary_inc_max,
                              "IncubateSecondary" = params$secondary_inc_max),
                            INCUBATION_STEPS )

  model <-
    model %>%
    add_constraint( inc_dur[step] >= inc_dur_min[step],
                    step = INCUBATION_STEPS) %>%
    add_constraint( inc_dur[step] <= inc_dur_max[step],
                    step = INCUBATION_STEPS)
  
  
  # ----- Exposure constraints -----

  exposure_steps <- new_subset( c("WashBeads",
                                  "WashSample",
                                  "WashPrimary",
                                  "WashSecondary",
                                  "AddEnhancer"),
                                PROCESSING_STEPS )
  
  post_exposure_steps <- new_subset( c("AddSample",
                                       "AddPrimary",
                                       "AddSecondary",
                                       "AddEnhancer",
                                       "Read"),
                                     PROCESSING_STEPS)
  
  exposure_bounds <- map_sets(exposure_steps, post_exposure_steps,
                           setNames( names(post_exposure_steps), names(exposure_steps) ) )
  
  exposures <- new_param( c("WashBeads" = params$wash_exposure_max,
                            "WashSample" = params$wash_exposure_max,
                            "WashPrimary" = params$wash_exposure_max,
                            "WashSecondary" = params$wash_exposure_max,
                            "AddEnhancer" = params$enhancer_exposure_max),
                          exposure_steps )
  
  model <-
    model %>%
    add_constraint( begin_time[plate, step2] - end_time[plate, step1] <= exposures[step1],
                    plate = PLATES,
                    step1 = exposure_steps,
                    step2 = post_exposure_steps,
                    # ompr removes all attributes, so must use a vectorized version
                    step2 == exposure_bounds(step1) )
  

  # ----- Resource conflict constraints -----

  # Assume that each physical resource (whether the pipettor/liquid handler,
  # washer or reader) can only process one plate at a time.


  # No conflicts are possible with only 1 plate...
  if (num_plates > 1L) {
    
    # Big M factor for enforcing the OR constraints
    bigM <- 1e6
    
    model <-
      model %>%
      
      ## Create the binary variables needed to enforce the OR constraint.
      add_variable( resOR[plate1, step1, plate2, step2],
                    plate1 = PLATES, step1 = PROCESSING_STEPS,
                    plate2 = PLATES, step2 = PROCESSING_STEPS,
                    # There's symmetry to the non-conflict constraints. This
                    # prevents the same constraints from being defined twice, and
                    # prevents a plate from being compared against itself:
                    plate1 < plate2,
                    # There's only a potential for conflict if the resources are
                    # the same (non-conflicting instruments should be able to
                    # run in parallel):
                    resources_map(step1) == resources_map(step2),
                    type = "binary" ) %>%
    
      # There's no resource conflict if the steps utilizing that resource don't
      # overlap in time. That is, there's no conflict if one step ends before the
      # other begins, OR if one begins after the other ends.
      #
      # The Big M notation ensures that if one constraint is violated, then the
      # other constrant MUST hold:
      add_constraint( end_time[plate1, step1] <= begin_time[plate2, step2] +
                        bigM * resOR[plate1, step1, plate2, step2],
                      plate1 = PLATES, step1 = PROCESSING_STEPS,
                      plate2 = PLATES, step2 = PROCESSING_STEPS,
                      plate1 < plate2,
                      resources_map(step1) == resources_map(step2)
                    ) %>%

      add_constraint( end_time[plate2, step2] <= begin_time[plate1, step1] +
                        bigM * (1 - resOR[plate1, step1, plate2, step2]),
                      plate1 = PLATES, step1 = PROCESSING_STEPS,
                      plate2 = PLATES, step2 = PROCESSING_STEPS,
                      plate1 < plate2,
                      resources_map(step1) == resources_map(step2)
                    )
  }
    
  
  # ----- Assay and run definitions -----
  
  # Order the assays by plate number:
  model <-
    model %>%
    add_constraint( assay_begin[plate1] <= assay_begin[plate2],
                    plate1 = PLATES, plate2 = PLATES,
                    plate1 < plate2 )
  
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
  
  # A run ends when all the assays have completed
  
  model <-
    model %>%
    add_constraint( assay_end[plate] <= run_end,
                   plate = PLATES )


  # ===== Objective =====
  
  # Minimize the overall run time for this number of plates
  model <- set_objective(model, run_end, sense = "min")

  
  return(model)
}


solve_schedule <- function(model) {
  
  solution <- solve_model(model, with_ROI(solver = "glpk"))
  
  return(solution)
}


repackage_solution <- function(model, solution, nplates = 1L) {
  
  PLATES <- seq_len(nplates)
  
  # Extract the solution
  
  ## Create the list of variables
  vars_list <- vars( run_end = NULL,
                     inc_dur = sets(Step = INCUBATION_STEPS),
                     assay_begin = sets(Plate = PLATES),
                     assay_end = sets(Plate = PLATES),
                     assay_dur = sets(Plate = PLATES),
                     begin_time = sets(Plate = PLATES, Step = STEPS),
                     end_time = sets(Plate = PLATES, Step = STEPS) )
  
  sols <- extract_solution(model, solution, vars_list)
  
  # Repackage this to be more user-friendly
  
  ## Incubation durations
  inc_dur <- rename(sols$inc_dur, Duration = value)
  
  ## Assay timings
  assay_begin <- rename(sols$assay_begin, StartTime = value)
  assay_end <- rename(sols$assay_end, EndTime = value)
  assay_dur <- rename(sols$assay_dur, Duration = value)
  
  assay_times <-
    inner_join(assay_begin, assay_end, by = "Plate") %>%
    inner_join(assay_dur, by = "Plate") %>%
    select(Plate, StartTime, EndTime, Duration) %>%
    arrange(Plate)
  
  ## Step timings
  step_begin <- rename(sols$begin_time, StartTime = value)
  step_end <- rename(sols$end_time, EndTime = value)
  
  step_times <-
    inner_join(step_begin, step_end, by = c("Plate", "Step")) %>%
    mutate( Duration = EndTime - StartTime) %>%
    select(Plate, Step, StartTime, EndTime, Duration) %>%
    arrange(Plate, Step)
  
  return( list(RunDuration = unname(sols$run_end), 
               IncubationDurations = inc_dur,
               AssayTimes = assay_times,
               StepTimes = step_times) )
}
