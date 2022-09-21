# load libraries

suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
})

base_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series"

minutes_since_last_update = function(file_name) {
  (as.numeric(as.POSIXlt(Sys.time())) -  
     as.numeric(file.info(file_name)$ctime)) / 60
}

load_data = function(file_name, column_name) {
  if(!file.exists(file_name) || 
     minutes_since_last_update(file_name) > 10) {
    data = read.csv(file.path(base_url, file_name), 
                    check.names = FALSE, stringsAsFactors = FALSE) %>%
      # select(-Lat, -Long) %>% 
      pivot_longer(-(1:2), names_to = "date", values_to = column_name)%>%
      mutate(
        date = as.Date(date, format = "%m/%d/%y"),
        `Province/State`=
          if_else(`Province/State` == "", "<all>", `Province/State`)
      )
    save(data, file = file_name)  
  } else {
    load(file = file_name)
  }
  return(data)
}