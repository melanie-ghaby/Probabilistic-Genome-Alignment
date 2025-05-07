#name: Mélanie Ghaby

import data_preprocessing
import query_generator_and_tester
import subarray_indentifier
import targeted_reprocessing
import diagonal_dictionary_recomposer
import alignment
import numpy as np
import time
import Levenshtein
import pandas as pd

#functions to fill out a matrix later used to determine the precision and
#correctness of the algorithm
def correct_subarray_identification(query_start_index, suggested_start, suggested_end, query):
    if query_start_index < suggested_start:
        return 0
    if query_start_index > suggested_end:
        return 0
    if (query_start_index + len(query)) > suggested_end:
        return 0
    return 1       

def correct_diagonal_identification(query_start_index, diagonal_start, diagonal_end, query):
    if query_start_index < diagonal_start:
        return 0
    if query_start_index > diagonal_end:
        return 0
    if (query_start_index + len(query)) > diagonal_end:
        return 0
    return 1      

def correct_start_index_identification(aligned_query, query_true_start, alignment_start):
    i = 0
    while i < len(aligned_query):
        if aligned_query[i] == '-':
            alignment_start += 1
            i += 1
        else: 
            i = 99999
            break
    
    return abs(alignment_start - query_true_start)
 

def correct_alignment_tester(aligned_query, generated_query):
    return (Levenshtein.ratio(aligned_query, generated_query) * 100)

#writing the computed information into an excel file for later analysis
def stats(time_matrix, correctness_matrix, word_sizes, file_name="tester1_results.xlsx"):
    time_col = ["Word Size", "Query Length", "Subarray Time", 
                    "Diagonal Range Time", "Alignment Time", "Total Time"]
    correctness_col = ["Word Size", "Query Length", "Subarray Correctness", 
                           "Diagonal Range Correctness", "Start Index Distance", 
                           "Levenshtein Ratio"]
    
    time_matrixdf = pd.DataFrame(time_matrix, columns=time_col)
    correctness_matrixdf = pd.DataFrame(correctness_matrix, columns=correctness_col)
    #writing to excel to later compute stats with R
    with pd.ExcelWriter(file_name, engine='openpyxl') as writer:
        time_matrixdf.to_excel(writer, sheet_name="Time Matrix", index=False)
        correctness_matrixdf.to_excel(writer, sheet_name="Correctness Matrix", index=False)


if __name__ == "__main__":
    num_iterations = 10
    time_matrix = np.zeros((num_iterations, 6)) 
    correctness_matrix = np.zeros((num_iterations, 6))  
    np.random.seed(100)

    i = 0
    correct_count = 0
    while i < num_iterations:
        word_size = np.random.randint(6, 11)  #randomly chosing a word size

        #data preprocessing
        nucleotide_array, probs_array = data_preprocessing.data_preprocessor(word_size)
        
        #Making the query to test the algorithm on
        query, real_start_index, intended_aligned_query = query_generator_and_tester.query_sequence_maker(nucleotide_array, probs_array)
        
        start_time = time.perf_counter()
        #identifying the initial subarray to explore
        #making a probabilistic dictionary for this subarray
        sub_start, sub_end = subarray_indentifier.first_pass_single_hits(nucleotide_array, query, word_size)
        prob_dict_filename = targeted_reprocessing.build_probabilistic_blast_database(nucleotide_array[sub_start:sub_end], probs_array[sub_start:sub_end], word_size, 0.95)
        step2_3time = time.perf_counter()

        #finding a slightly longer than the query range to align the query
        #CAREFUL: INDECES NOW HAVE OFFSET DUE TO SELECTED RANGE
        start_alignment, end_alignment = diagonal_dictionary_recomposer.diagonal_dictionary_maker_cleaner(word_size, query, prob_dict_filename, sub_start)
        step4_time = time.perf_counter()
        
        #aligning the query sequence with the region of the database
        opt_score, aligned_query, aligned_dtb = alignment.full_algorithm(query, nucleotide_array[start_alignment: end_alignment], probs_array[start_alignment: end_alignment])
        end_time = time.perf_counter()
  
        time_matrix[i, 0] = word_size  
        time_matrix[i, 1] = len(query)
        time_matrix[i,2] = step2_3time-start_time
        time_matrix[i,3] = step4_time-step2_3time
        time_matrix[i,4] = end_time - step4_time
        time_matrix[i,5] = end_time - start_time

    
        correctness_matrix[i, 0] = word_size
        correctness_matrix[i, 1] = len(query)
        correctness_matrix[i, 2] = correct_subarray_identification(real_start_index, sub_start, sub_end, query)        
        correctness_matrix[i, 3] = correct_diagonal_identification(real_start_index, start_alignment, end_alignment, query)
        correctness_matrix[i, 4] = correct_start_index_identification(aligned_query, real_start_index, start_alignment)
        #stripping them both of gaps in front and at the end since aligned query has inserted gaps
        #since aligned over its own size
        correctness_matrix[i, 5] = correct_alignment_tester(aligned_query.strip('-'), intended_aligned_query.strip('-'))
        print(i)
        i += 1
        

    stats(time_matrix, correctness_matrix, "tester1_results.xlsx")
    
    

