source("scheduler.R")
source("gantt-chart.R")


test_model <- function(nplates = 1L, by = "plate") {
  
  model <- build_model(num_plates = nplates, regular = TRUE)
  solution <- solve_schedule(model)
  sol <- repackage_solution(model, solution, nplates)
  
  # resources <-  data.frame( 
  #   Step = c("AddBeads",
  #            "WashBeads", 
  #            "AddSample",
  #            "IncubateSample",
  #            "WashSample",
  #            "AddPrimary",
  #            "IncubatePrimary",
  #            "WashPrimary",
  #            "AddSecondary",
  #            "IncubateSecondary",
  #            "WashSecondary",
  #            "AddEnhancer",
  #            "Read"),
  #   Resource = factor(
  #     c("Pipettor",
  #       "Washer",
  #       "Pipettor",
  #       "Incubation",
  #       "Washer",
  #       "Pipettor",
  #       "Incubation",
  #       "Washer",
  #       "Pipettor",
  #       "Incubation",
  #       "Washer",
  #       "Pipettor",
  #       "Reader"),
  #     levels = c("Pipettor", "Washer", "Reader", "Incubation")),
  #   stringsAsFactors = FALSE
  # )
  
  # timings <- 
  #   sol$StepTimes %>%
  #   # Convert "Step" to a character vector, to keep ggplot2 from complaining,
  #   # and to ensure that the join aligns properly:
  #   mutate(Step = as.character(Step)) %>%
  #   rename(Start = StartTime, Stop = EndTime) %>%
  #   inner_join(resources, by = "Step")
  
  resources <- new_set(c("Pipettor", "Washer", "Reader", "Incubation"))
  
  # Map resources 
  res_map <- map_sets( STEPS, resources, 
                             c(  "AddBeads"          = "Pipettor",
                                 "WashBeads"         = "Washer",
                                 "AddSample"         = "Pipettor",
                                 "IncubateSample"    = "Incubation",
                                 "WashSample"        = "Washer",
                                 "AddPrimary"        = "Pipettor",
                                 "IncubatePrimary"   = "Incubation",
                                 "WashPrimary"       = "Washer",
                                 "AddSecondary"      = "Pipettor",
                                 "IncubateSecondary" = "Incubation",
                                 "WashSecondary"     = "Washer",
                                 "AddEnhancer"       = "Pipettor",
                                 "Read"              = "Reader"
                             )
  )
  
  timings <-
    sol$StepTimes %>%
    mutate( Resource = factor( names(res_map(as.character(Step))),
                               levels = names(resources) ) ) %>%
    filter( Resource != "Incubation") %>%
    rename(Start = StartTime, Stop = EndTime)
  
  show(gantt_chart(timings, by = by))
  
  return(sol)
}
