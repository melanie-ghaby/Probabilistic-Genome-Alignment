#name: Mélanie Ghaby


############################### PART 2 ########################################
#Divides the entire database into subarrays, and decides to proceed
#with the algorithm on only one subarray depending on the one that contains the
#greatest number of hits

import json
import math
import matplotlib.pyplot as plt

#returns a list of tuple of hits in the following form; 
#(starting position of hit in dtb, starting position of hit in query)
#this list is sorted on ascending order of hits in database
def find_query_word_hits(query, word_size, database_filename):
    with open(database_filename, "r") as file:
        database = json.load(file)
    
    positions_ofhits = []
    #extracting words from the query sequence
    for i in range(len(query) - word_size + 1):
        word = query[i:i + word_size]  
        #for word in the database, pos is the position in database, i in query
        '''
        if word in database:
            positions_ofhits += [(current_pos, i) for current_pos in database.get(word)]
        '''
        if word in database:
            cur_positions = database.get(word)
            for current_pos in cur_positions:
                # Add the tuple (current_pos, i) to the 'positions_ofhits' list
                positions_ofhits.append((current_pos, i))
        
    #sort list based based on ascending order of database positions 
    positions_ofhits.sort(key=lambda x:x[0])
    return positions_ofhits



#helper function later used to identify the subarray of the database
# that has the most common words with the query
#given a list, it counts the number of occurences
#of a value between two values
def count_occurrences_in_range(pos_list, start, finish):
    count = 0
    query_word_hits = []
    for index_dtb, index_query in pos_list:  
        if start <= index_dtb< finish:
            count += 1
            query_word_hits.append(index_query)
        if index_dtb >= finish:
            break  #pos_list sorted in ascending order of index_dtb, so can break
            
    #number of unique query hits within the range [start, finish]
    unique_indices_hits = len(set(query_word_hits))

    return (count* unique_indices_hits)


#returns tuples in the following format 
#(start_pos in database, end_post in database, #k_common_words)
#divides the list of positions of common words between query sequence and 
#database into chunks of equal sizes 
def subarrays_tuples_positions(pos_list,nucleotide_array, query):
    subarrays = math.floor(len(nucleotide_array)/(len(query)*5))
    first_delim = len(nucleotide_array) /subarrays  
    tuples = []

    #counting occurences for subarrays 
    for i in range(subarrays):
        start = math.floor(i * first_delim)
        end = math.ceil((i + 1) * first_delim)
        if (end > len(nucleotide_array)):
            end = len(nucleotide_array) - 1
        occurrences = count_occurrences_in_range(pos_list, start, end)
        tuples.append((start, end, occurrences))
    
    return tuples

#finds max subbarray hit count, and joins it with left most and right most
#if they exidt
def max_common_words_finder(tuples, query):
    #reverse = True for descending order to find the max
    sorted_tuples = sorted(tuples, key=lambda x: x[2], reverse=True)
    max_tuple = sorted_tuples[0]
    max_index = tuples.index(max_tuple)
    left_tuple = tuples[max_index - 1] if max_index > 0 else None
    right_tuple = tuples[max_index + 1] if max_index < (len(tuples) - 1) else None
    start_index = left_tuple[0] if left_tuple else max_tuple[0]
    end_index = right_tuple[1] if right_tuple else max_tuple[1]
    total_occurrences = (left_tuple[2] if left_tuple else 0) + max_tuple[2] + (right_tuple[2] if right_tuple else 0)
    first_subarray = (start_index, end_index, total_occurrences)
    
    return first_subarray

#max subarray and its left most and right most are redivded into smaller ones
#and heir occurences count are recounted
def subarrays_refiner(pos_list,subarray_range, query):
    current_subarray = subarray_range[1]- subarray_range[0]
    offset=subarray_range[0]
    subarrays = math.floor(current_subarray/(len(query)*2.5))
    first_delim = (current_subarray / subarrays) 
    tuples = []
    
    for i in range(subarrays):
        start = math.floor(i * first_delim) + offset
        end = math.ceil((i + 1) * first_delim) + offset
        if (end > subarray_range[1]):
            end = subarray_range[1] 
        occurrences = count_occurrences_in_range(pos_list, start, end)
        tuples.append((start, end, occurrences))
    
    return tuples

    


#returns starting and ending positions of chunks in the database where the
#query sequence is most likely to originate from, inspired by the most commun 
#number of words

def first_pass_single_hits(nucleotide_array, query, word_size):
    positions = find_query_word_hits(query, 8, 'w_8.json')
    #visualize_subarrays_with_counts(positions, nucleotide_array, query)
    tuples = subarrays_tuples_positions(positions, nucleotide_array, query)
    highest =  max_common_words_finder(tuples, query)
    new_tuples = subarrays_refiner(positions, highest, query)
    new_highest = max_common_words_finder(new_tuples, query)
    return (int(new_highest[0]), (int(new_highest[1])+word_size) )


###################subarray visualizer#######################
#used for the report part
def visualize_subarrays_with_counts(pos_list, nucleotide_array, query):
    tuples = subarrays_tuples_positions(pos_list, nucleotide_array, query)
    subarray_numbers = list(range(len(tuples)))  
    occurrences = [t[2] for t in tuples]  
    
    plt.figure(figsize=(10, 6))
    
    plt.plot(subarray_numbers, occurrences, marker='o', color='black', linestyle='-', markersize=6)
    
    plt.title('Rewarded Unique Hits Based on Subarray Number')
    plt.xlabel('Subarray Number')
    plt.ylabel('Number of Hits * Uniquess Reward')
    
    plt.grid(True)
    plt.show()

    


            
        

