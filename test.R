source("scheduler.R")
source("gantt-chart.R")


test_model <- function(nplates = 1L) {
  
  model <- build_model(num_plates = nplates)
  solution <- solve_schedule(model)
  sol <- repackage_solution(model, solution, nplates)
  
  resources <-  data.frame( 
    Step = c("AddBeads",
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
             "Read"),
    Resource = factor(
      c("Pipettor",
        "Washer",
        "Pipettor",
        "Incubation",
        "Washer",
        "Pipettor",
        "Incubation",
        "Washer",
        "Pipettor",
        "Incubation",
        "Washer",
        "Pipettor",
        "Reader"),
      levels = c("Pipettor", "Washer", "Reader", "Incubation")),
    stringsAsFactors = FALSE
  )
  
  timings <- 
    sol$StepTimes %>%
    # Convert "Step" to a character vector, to keep ggplot2 from complaining,
    # and to ensure that the join aligns properly:
    mutate(Step = as.character(Step)) %>%
    rename(Start = StartTime, Stop = EndTime) %>%
    inner_join(resources, by = "Step")
  
  gantt_chart(timings, by = "plate")
}
