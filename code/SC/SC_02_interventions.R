# identifying dates which slope can be calculated 
# # this data will be used for filtering the samples before calculating 
# # the slope to avoid inclusion of any intervention that affects the slope
# # (Tolvaptan usage, nephrectomy, kidney transplant, dialysis)

mct_table$HD[which(mct_table$HD == "continuos gene syndrome")] = NA
problems_mct = unique(mct_table[unique(c(which(!is.na(mct_table$HD)), 
                                              grep("Nephrekt|ntx|einzel",mct_table$Comment, 
                                                   ignore.case = T))),"identifier"])

require(stringr)
commtfix = lapply(str_split(mct_table[problems_mct,"Comment"], " "), function(x) gsub(",", "", x))

mct_table$dates = NA
mct_table[problems_mct,"dates"] = unlist(lapply(commtfix, function(x) grep("\\d|Einzel", x, value = T)[1]))
mct_table$dates[grep("TKV|Einzelniere", mct_table$Comment)] = "01.01.1900"
mct_table$dates[which(nchar(mct_table$dates) < 5)] = paste("01.01.", mct_table$dates[which(nchar(mct_table$dates) < 5)], sep = "")

require(lubridate)
mct_table$dates = as.Date(parse_date_time(mct_table$dates, orders = c("dmy")))

mct_table$HDdate = gsub(".* ", "", mct_table$HD)
setdate1900 = intersect(which(!is.na(mct_table$HDdate)), 
                        grep("\\d", (gsub(".* ", "", mct_table$HDdate)), 
                             invert = T))
mct_table$HDdate[setdate1900] = "01.01.1900"

require(janitor)
mct_table$coeldates = as.Date(NA)
mct_table$coeldates[nchar(mct_table$HDdate) > 5 & !is.na(mct_table$HDdate)] = as.Date(mct_table$HDdate[nchar(mct_table$HDdate) > 5 & !is.na(mct_table$HDdate)], format = "%d.%m.%Y")
mct_table$coeldates[which(nchar(mct_table$HDdate) == 5)] = janitor::excel_numeric_to_date(as.numeric(mct_table$HDdate[which(nchar(mct_table$HDdate) == 5)]))
mct_table$coeldates = as.Date(mct_table$coeldates)

mct_table$DATES = as.Date(NA)
counter = 1
invisible(apply(mct_table[,c("dates","coeldates")],1, function(x){
  i = which.min(as.Date(x))
  mct_table$DATES[counter] <<- x[i]
  counter <<- counter + 1
}))

tolv_used_p = mct_table[which(mct_table$Tolvaptan == TRUE | !is.na(mct_table$DATES)),]
tolv_used_p$`Start 45/15` = as.Date(tolv_used_p$`Start 45/15`)

tolv_used_p$start4515orginal = tolv_used_p$`Start 45/15`
tolv_used_p$`Start 45/15` = as.Date(NA)
counter = 1
invisible(apply(tolv_used_p[,c("start4515orginal","DATES")], 1, function(x){
  i = which.min(as.Date(x))
  tolv_used_p$`Start 45/15`[counter] <<- x[i]
  counter <<- counter + 1
}))

tolv_used_p$patient_id = tolv_used_p$identifier
tolv_used_p2 = tolv_used_p
