
### Find the JSON path information in the appropriate directory.
jinfo <- "COAD_ra.json"
if (!file.exists(jinfo)) stop("Cannot locate file: '", jinfo, "'.\n", sep='')
### parse it
library(rjson)
temp <- fromJSON(file = jinfo)
paths <- temp$paths
detach("package:rjson")
### clean up
rm(jinfo, temp)