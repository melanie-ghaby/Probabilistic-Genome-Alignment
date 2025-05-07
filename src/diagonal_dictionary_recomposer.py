#name: Mélanie Ghaby


import json
import statistics



################################### PART 4 ###################################
#Identifies the range of diagonals for which
#the number of hits is great enough to be the potential originator of the query


#makes a list of tuples of the form:
#(starting position in database,  starting pos in query, probability of word correctness)
#the start_index_offset corresponds to the starting_index in the database of the subarray
#chosen to build the dictionary
def find_query_word_hits_with_probs(query, word_size, prob_database_file, start_index_offset):
    with open(prob_database_file, "r") as file:
        database = json.load(file)
    
    positions = []
    #extracting words from the query sequence
    for i in range(len(query) - word_size + 1):
        word = query[i:i + word_size]  
        
        #for every word found in the database
        #add tuples (pos in dtb, pos i in query , prob of word, word itself) to the list
        #take into account the offset of taking a subarray
        #using list comprehenion to makes things quicker
        if word in database:
            positions += [(pos[0] + start_index_offset, i, pos[1], word) for pos in database.get(word)]
    
    #sort list based on ascending database positions 
    positions.sort(key=lambda x: x[0])
    return positions


#makes a dictionary where a key is a diagonal and its values are the tuples of the form
# (position in database, position in query, probs of word correctness, word, word_size)
def diagonal_dictionary_maker(positions_list, word_size):
    diagonal_dict = {} 
    #computing the diagonal
    for db_pos, query_pos, prob, word in positions_list:
        diagonal = db_pos - query_pos  
        
        #adding every tuple to its respective diagonal 
        if diagonal not in diagonal_dict:
            diagonal_dict[diagonal] = []  
        diagonal_dict[diagonal].append((db_pos, query_pos, prob, word, word_size))
    return diagonal_dict


#these are two helper functions 
#using list comprehension requires less function calls
def remove_single_hit_diagonals(diago_dict):
    filtered_diago_dict = {diagonal:  hits for diagonal, hits in diago_dict.items() if len(hits) >= 2}
    #used dictionary comprehension for making things fasetr
    return filtered_diago_dict

def order_diagonal_keys(diago_dict):
    sorted_diago_dict = {k: diago_dict[k] for k in sorted(diago_dict.keys())}
    return sorted_diago_dict


#returns the metrics of the dictionary of diagonals in a dictionary format
def diagonal_metrics(diago_dict):
    diagonal_metrics_dict = {}
    #going through every diagonals and its hits
    for diagonal, words_list in diago_dict.items():
        #diagonals with less than 2 hits should be filtered out by now
        #computing the average distance between hits next to each other
        #list comprehension to speed up process
        positions = [word_tuple[1] for word_tuple in words_list]
        distances = [positions[i + 1] - positions[i] for i in range(len(positions) -1)]
        avg_distance = sum(distances) /len(distances)
        numb_hits = len(words_list)

        diagonal_metrics_dict[diagonal] = {'avg_distance': avg_distance, 'numb_hits': numb_hits}

    return diagonal_metrics_dict


#returns a dictionary containing the diagonals with the highest number of hits only.
#the number of diaognals returning is a parameter that is 15 by default
#if there is a tie at the last position, adds them all
def select_top_diagonals_by_hits(diago_dict, top_n=15):
    #computing metrics
    diagonal_metrics_dict = diagonal_metrics(diago_dict)
    
    #sorting dictionary by ascending number of hits per diagonal for the metrics dictionary
    sorted_diagonals = sorted(diagonal_metrics_dict.items(), key=lambda x: x[1]['numb_hits'], reverse=True)
    
    #finding cutoff to make it to the top_n diagonals
    top_diagonals = {}
    if len(sorted_diagonals) >= top_n:
        last_hit_count = sorted_diagonals[top_n - 1][1]['numb_hits']
    else:
        last_hit_count = sorted_diagonals[-1][1]['numb_hits']

    #add diagonals with more hits than the cutoff 
    for diagonal, metrics in sorted_diagonals:
        if metrics['numb_hits'] >= last_hit_count:
            top_diagonals[diagonal] = diago_dict[diagonal]
     
    #reconsider bounding numbre later
    '''
    #if tying for last position, pick the one closest to the max hits diagonal
    #if there are more than allowed
    if len(top_diagonals) > top_n:
        #diagonal with the maximum number of hits found using the key and the max function
        max_hits_diagonnal = max(diagonal_metrics_dict, key=lambda k: diagonal_metrics_dict[k]['numb_hits'])
        top_diago_list = list(top_diagonals.items())
        #sorting the diagonals by proximity to the one with max hits
        top_diago_list.sort(key=lambda x:abs(x[0] -max_hits_diagonnal))
        #selecting top n diagonals from the sorted one based on proximity
        top_diago = dict(top_diago_list[:top_n])
    '''  
    
    return top_diagonals



#computes the average distance between diagonals, and removes those whose 
#distance to the next diagonal is greater than the average distance
#between consecutive diagonals or median 
#also removes diagonals that are more than 15 indexes away from the diagonal
#with the maximum number of hits
#the diagonal with the maximum number of hits is always included
def filter_diagonals_by_avg_distance(diago_dict):
    #comuting diagonal metrices
    diagonal_metrics_dict = diagonal_metrics(diago_dict)
    #finding key corresponding to the diagonal with the max number of hits
    #this is in numb hits
    max_hits_diagonal = max(diagonal_metrics_dict, key=lambda k: diagonal_metrics_dict[k]['numb_hits'])
    
    #sorting all diagonal in ascending order (so smalles to biggest diago)
    diagonals = sorted(diago_dict.keys())
    #computing distances between pairs of diagonals
    #using list comprehnion for more efficient computations
    distances = [diagonals[i + 1] - diagonals[i] for i in range(len(diagonals) - 1)]
    
    
    #computing average distance and median distance
    avg_distance = 0
    median_distance=0
    if distances:
        avg_distance = sum(distances) / len(distances) 
        median_distance = statistics.median(distances) 

    filtered_diagonals = []
    #including the diagonal with the max number of hits always
    filtered_diagonals.append(max_hits_diagonal)
    
    for i in range(1, len(diagonals)):
        diagonal = diagonals[i]
        
        #only keep diagonals within a median and avg distance from one another 
        #in order to prevent outliers
        if (diagonal - diagonals[i - 1] <= median_distance  and 
            diagonal - diagonals[i - 1] <= avg_distance and 
            abs(diagonal - max_hits_diagonal) <= 15):
            filtered_diagonals.append(diagonal)
    #recreate new dictionary but this time only with top diagonlas
    filtered_diago_dict = {diagonal: diago_dict[diagonal] for diagonal in filtered_diagonals}

    return filtered_diago_dict

#returns the first and last key of entries of a dictionary
def first_and_last_key(diago_dict):
    keys_list = list(diago_dict.keys())
    return (keys_list[0], keys_list[-1])


def diagonal_dictionary_maker_cleaner(word_size, query, prob_database_filename,start_index_offset):
    positions = find_query_word_hits_with_probs(query, word_size, prob_database_filename, start_index_offset)
    diagonal_dictionary = diagonal_dictionary_maker(positions,word_size)
    without_one = remove_single_hit_diagonals(diagonal_dictionary)
    ordered_dictio = order_diagonal_keys(without_one)
    relevant_diagonals = select_top_diagonals_by_hits(ordered_dictio)
    order_relevant = order_diagonal_keys(relevant_diagonals)
    remove_too_far_diago = filter_diagonals_by_avg_distance(order_relevant)
    order_relevant2 = order_diagonal_keys(remove_too_far_diago)
    (start,end) = first_and_last_key(order_relevant2)
    return (int(start)-10, (int(end) + len(query)+10))

