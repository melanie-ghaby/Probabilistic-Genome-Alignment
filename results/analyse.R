#Name: Melanie Ghaby

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

official_dataset <- "C:/Users/ghaby/OneDrive - McGill University/Desktop/Fall 2024/COMP 561/official_dataset.xlsx"

##########The data here is not well formatted, this is done in the final report manually##########

##############Time analysis###################
data <- read_excel(official_dataset, sheet = "Time Matrix")

#general summary statistics
time_gen_stats <- data %>%
  summarise(
    Avg_subarray_processing_time = mean(`Subarray Time`),
    median_subarray_processing_time = median(`Subarray Time`),
    Avg_diagonal_range_finder_time = mean(`Diagonal Range Time`),
    Median_diagonal_range_finder_time = median(`Diagonal Range Time`),
    Avg_alignment_time = mean(`Alignment Time`),
    Median_alignment_Time = median(`Alignment Time`),
    Avg_full_time = mean(`Total Time`),
    Median_full_time = median(`Total Time`),
  )

summary_stats_long <- time_gen_stats %>%
  pivot_longer(cols = everything(), names_to = "Corresponding statistic", values_to = "Value")

summary_stats_long


#make intervals for query lengths
data <- data %>%
  mutate(query_lengths = cut(`Query Length`,
                                breaks = c(0, 500, 1000, 1500, 2000, 2510),
                                labels = c("248-500", "500-1000", "1000-1500", "1500-2000", "2000-2510"),
                                include.lowest = TRUE))


#visualizing the trend to see total running time vs query length and word size
ggplot(data, aes(x = `Query Length`, y = `Total Time`, group = factor(`Word Size`))) +
  geom_point(size = 2, alpha = 0.7) + facet_wrap(~`Word Size`, scales = "free_y") + 
  labs(
    title = "Total Running Time vs. Query Length",
    x = "Query Length (nucleotides)",
    y = "Total Running Time (s)"
  ) + theme_minimal() 



##################correctness analysis############

data <- read_excel(official_dataset, sheet = "Correctness Matrix")

#that is because if an initial step is incorrect, the rest is certainly incorrect
data[ data$`Diagonal Range Correctness` == 0, 5:ncol(data) ] <- NA
data[ data$`Subarray Correctness` == 0, 4:ncol(data) ] <- NA

#make intervals for query lengths
data <- data %>%
  mutate(query_lengths = cut(`Query Length`,
                             breaks = c(0, 500, 1000, 1500, 2000, 2510),
                             labels = c("248-500", "500-1000", "1000-1500", "1500-2000", "2000-2510"),
                             include.lowest = TRUE))

#checking the total number of simulations performed for each category
numbers <- data %>%
  group_by(`Word Size`, query_lengths) %>%
  summarise(
     number_conducted =n()
  ) 

numbers 

#some statistics about the success of the subarray identification
# and the diagonal range identification
result_1 <- data %>%
  group_by(`Word Size`, query_lengths) %>%
  summarise(
    subarray_success = sum(`Subarray Correctness` == 1, na.rm = TRUE),
    diagonal_success = sum(`Diagonal Range Correctness` == 1, na.rm = TRUE),
    total_subarray = sum(!is.na(`Subarray Correctness`)),  
    total_diagonal = sum(!is.na(`Diagonal Range Correctness`)),  
    subarray_fraction = paste(subarray_success, total_subarray, sep = "/"),  
    diagonal_fraction = paste(diagonal_success, total_diagonal, sep = "/"),  
    subarray_ration = (subarray_success/total_subarray),
    diagonal_ration = (diagonal_success/total_diagonal)
  )  %>%
  select(`Word Size`, query_lengths,  subarray_ration,  diagonal_ration)

print(result_1)

#stats about the distance metrics from real query start
result_2 <- data %>%
  group_by(`Word Size`, query_lengths) %>%
  summarise(
    avg_start_index_distance = mean(`Start Index Distance`, na.rm = TRUE),
    num_start_index_distance_perfect = sum(`Start Index Distance` == 0, na.rm = TRUE),
    total_start_index_distance_tests = sum(!is.na(`Start Index Distance`)),
    perfect_ration =( num_start_index_distance_perfect /total_start_index_distance_tests)*100,
    worst_ditance = max(`Start Index Distance`, na.rm = TRUE),
  ) %>% select(c(`Word Size`, query_lengths, "avg_start_index_distance", "perfect_ration", worst_ditance))

print(result_2)

#Since there seems to be no pattern in relation of word size of query length, the
#distance from the query start will be average over all data points: 
data %>% summarise(mean_start_index_distance = mean(`Start Index Distance`, na.rm = TRUE))


#levenshtein ration based on query lengths
ggplot(data, aes(x = `Query Length`, y = `Levenshtein Ratio`, col = as.factor(`Word Size`))) + 
  geom_point(na.rm = TRUE) +
  labs(title = "Levenshtein Ratio vs. Query Length",
       x = "Query Length (nucleotides)",
       y = "Levenshtein Ratio (%)",
       color = "Word Size") +
  theme_minimal()

