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
    Resource = c("Pipettor",
                 "Washer",
                 "Pipettor",
                 "Incubate",
                 "Washer",
                 "Pipettor",
                 "Incubate",
                 "Washer",
                 "Pipettor",
                 "Incubate",
                 "Washer",
                 "Pipettor",
                 "Reader"),
    stringsAsFactors = FALSE
  )
  
  timings <- 
    sol$StepTimes %>%
    rename(Start = StartTime, Stop = EndTime) %>%
    inner_join(resources, by = "Step")
  
  gantt_chart(timings, by = "plate")
}
