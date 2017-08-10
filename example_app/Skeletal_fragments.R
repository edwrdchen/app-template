#read data

nlp_data <- read.table("sentences_nlp352.txt", fill=TRUE, na.strings=c("", "NA"))

rock_units <- read.csv("https://macrostrat.org/api/units?all&format=csv")

rock_units$unit_name <- gsub("Mbr", "", rock_units$unit_name)
rock_units$unit_name <- gsub("Member", "", rock_units$unit_name)
rock_units$unit_name <- gsub("Fm", "", rock_units$unit_name)
rock_units$unit_name <- gsub("Gp", "", rock_units$unit_name)

time_scale <- read.csv("https://paleobiodb.org/data1.2/intervals/list.txt?all_records")

#Define functions
number_extract <- function(dependencies) {
  
  subset(dependencies["dependent"], dependencies["type"] == "nummod" & dependencies["governor"] == "%")
  
}


#____________________________________________________________________________________________________________

formation_name_extract <- function(dependencies) {
  
  
  formation_name_unprocessed <- c(subset(as.character(dependencies$dependent), dependencies$type=="compound"), subset(as.character(dependencies$governor), dependencies$type=="compound"))
  
  formation_name_samplefree <- gsub("Sample", "", formation_name_unprocessed)
  formation_name_samplesfree <- gsub("Samples", "", formation_name_samplefree)
  formation_name_withtime <- unique(formation_name_samplesfree[grep("[[:upper:]]", strtrim(formation_name_samplesfree, 1))])
  timePeriodNames <- time_scale$interval_name
  subset(formation_name_withtime, !formation_name_withtime %in% timePeriodNames)

}

#___________________________________________________________________________________________


time_period_extract <- function(dependencies) {
  
  formation_name_unprocessed <- c(subset(as.character(dependencies$dependent), dependencies$type=="compound"), subset(as.character(dependencies$governor), dependencies$type=="compound"))
  formation_name_samplefree <- gsub("Sample", "", formation_name_unprocessed)
  formation_name_samplesfree <- gsub("Samples", "", formation_name_samplefree)
  formation_name_withtime <- unique(formation_name_samplesfree[grep("[[:upper:]]", strtrim(formation_name_samplesfree, 1))])
  timePeriodNames <- time_scale$interval_name
  time_periods <- intersect(timePeriodNames, unlist(strsplit(paste(formation_name_withtime, collapse = " "), ' ')))
  if (length(time_periods) > 1) {return(NA)} else {time_periods}
   
}

#_____________________________________________________________________________________________________________________


time_extract <- function(formation_name) {

  formation_name <- gsub("[Mm]ember", "", formation_name)
  formation_name <- gsub("[Mm]br", "", formation_name)
  formation_name <- gsub("[Mb]b", "", formation_name)
  
  formation_name <- gsub("[Ff]ormation", "", formation_name)
  formation_name <- gsub("[Ff]m", "", formation_name)
  
  formation_name <- gsub("[Gg]roup", "", formation_name)
  formation_name <- gsub("Gp", "", formation_name)

  rock_time <- subset(rock_units$b_age, rock_units$unit_name == formation_name)

  median(rock_time, na.rm = TRUE)
}


NLP_Info_Extract <- function(nlp) {
  tryCatch({
    #Extracts dependent words from single sentence
    clean_commas <- gsub("\",\"", "|", nlp[[4]], fixed = TRUE)
    if (grepl("\"[0-9]+,[0-9]+\"", clean_commas)==TRUE) {
      sentence <- as.character(nlp[[4]])
      words.df <- read.csv(text=sentence, header=F)
      words <- apply(words.df, 1, function(x) gsub("[{}]", "", x))
    } else {
      cleanWord <- lapply(clean_commas, function(x) gsub("[{}]", "", x))
      words <- unlist(strsplit(as.character(cleanWord), "\\,"))
    }
    #Extract type relationship words
    nlp_parts <- strsplit(as.character(nlp[[8]]), "\\,")
    type <- unlist(lapply(nlp_parts, function(x) gsub("[{}]", "", x)))

    #Extract dependency numbers
    cleanNumber <- strsplit(as.character(nlp[[9]]), "\\,")
    dependency_number <- unlist(lapply(cleanNumber, function(x) gsub("[{}]", "", x)))
    
    dependency_number[which(dependency_number==0)] <- NA
    
    dependencies <- data.frame("dependent"=words, "type"=type)
    
    dependencies$governor <- dependencies$dependent[as.numeric(dependency_number)]
    
    dependencies
  },
  error=function(error_message){
    message("Incomplete final line found by readTableHeader on 'text'")
    message(error_message)
    return(NA)
    #Extracts dependent words from single sentence
    sentence <- as.character(nlp[[4]])
    words.df <- read.csv(text=sentence, header=F)
    words <- apply(words.df, 1, function(x) gsub("[{}]", "", x))
    
    #Extract type relationship words
    parts <- as.character(nlp[[8]])
    type.df <- read.csv(text=parts, header=F)
    type <- apply(type.df, 1, function(x) gsub("[{}]", "", x))
    
    #Extract dependency numbers
    numb <- as.character(nlp[[9]])
    number.df <- read.csv(text=numb, header=F)
    dependency_number <- apply(number.df, 1, function(x) gsub("[{}]", "", x))
    
    dependency_number[which(dependency_number==0)] <- NA
    
    dependencies <- data.frame("dependent"=words, "type"=type)
    
    dependencies$governor <- dependencies$dependent[as.numeric(dependency_number)]
    
    dependencies
    }
  )
}


#find sentences that contain keywords of interest
sentences_no <- unique(c(grep("[Ss]keletal", nlp_data$V4), grep("[Ss]hell", nlp_data$V4), grep("[Ff]ossil", nlp_data$V4)))

good_sentences_blank <- nlp_data[sentences_no,]

good_sentences <- good_sentences_blank[complete.cases(good_sentences_blank),]



#create data frame with governor, dependent, and relationship type
sentence_parse <- apply(good_sentences, 1, function(x) NLP_Info_Extract(x))

#extracts time periods from previous data frame
time_period <- lapply(sentence_parse, time_period_extract)

time_period[lengths(time_period)==0] <- NA #converts empty to NA

#extracts formation names
formations <- lapply(sentence_parse, formation_name_extract)

formation_clean <- sapply(formations, function(x) paste(lapply(x, paste, collapse=", "), collapse=" "))


#finds skeletal abundance
skeletal_abund <- number_extract(sentence_parse)

skeletal_abund <- subset(skeletal_abund, skeletal_abund$sentence_no %in% rownames(table(skeletal_abund$sentence_no))[which(table(skeletal_abund$sentence_no) == 1)])


#finds age of formation in millions of years
age_data <- sapply(formation_clean, time_extract)

final_results <- data.frame(formation_clean, time_period=unlist(time_period), age_data)

#adds skeletal grains 
final_results$skeletal <- skeletal_abund$number_only[match(rownames(final_results), skeletal_abund$sentence_no)]


subset(final_results, is.na(final_results$age_data) == F | is.na(final_results$time_period)==F)


