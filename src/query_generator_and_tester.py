#name: Mélanie Ghaby


#############################TESTING AREA #####################################
#This is where a query sequence is generated according to the original sequence
#And the respective probabilities of correctness of its nucleotides, and it is 
#Then tested on the algorithm in order to determin its correctness


import numpy as np

np.random.seed(100)
#A mutation rate of 0.1 and an indel rate of 0.05 were chosen as an approximation of
#the actual rates in real life that I do not know about

#This generates a random sequence according ot the probability of each nucleotide
#Being in a certain position. It also mutates some nucleotides and indels
def generate_query_sequence(nucleotide_array, probs_array, mutation_rate=0.1, indel_rate=0.05):
    #choosing a random query length between 150 and 2500 nucleotides
    query_length = np.random.randint(150,2500)
    #choosing a random starting index
    start_index = np.random.randint(0, len(nucleotide_array) - query_length)
    
    #generating a query according to the probabilities imported
    query_sequence = []
    for i in range(query_length):
        position = start_index + i
        nucleotide_probs = probs_array[position]
        most_likely_nucleotide = nucleotide_array[position]
        incorrect_prob = (1 - nucleotide_probs) / 3  
        adjusted_probs = np.array([
            nucleotide_probs if nuc == most_likely_nucleotide else incorrect_prob
            for nuc in ['A', 'C', 'G', 'T']])
        
        #making sure the probabilities sum to 1 incase of mismatch due to decimals
        adjusted_probs /= adjusted_probs.sum()
        #adding a nucleotide to the sequence according to the probability of each
        #at that position
        nucleotide_choice = np.random.choice(['A', 'C', 'G', 'T'], p=adjusted_probs)
        query_sequence.append(nucleotide_choice)
    
    query_sequence = ''.join(query_sequence)
    query_sequence = list(query_sequence)
    
    #mutating some nucleotides
    for i in range(len(query_sequence)):
        if np.random.rand() < mutation_rate:
            query_sequence[i] = np.random.choice(['A', 'C', 'G', 'T'])

    #indels
    final_query = []
    mapped_query = [] #this one will store the way the mapped query should look like
    i = 0
    while i < len(query_sequence):
        #case of indel happening
        if np.random.rand() < indel_rate:
            #1 is for insertion and 0 for deletion
            insert_or_delete = np.random.choice([0, 1])
            
            if insert_or_delete == 0:
                #skipping adding the nucleotide at the current position since deleted
                mapped_query.append('-');
                i += 1
                continue
            else:
                #inserting a random nucleotide to the sequence
                inserted_nuc = np.random.choice(['A', 'C', 'G', 'T'])
                final_query.append(inserted_nuc)
                mapped_query.append(inserted_nuc)
                continue
        final_query.append(query_sequence[i])
        mapped_query.append(query_sequence[i])
        i += 1
               
    return ''.join(final_query), start_index, ''.join(mapped_query)


def query_sequence_maker(nucleotide_array, probs_array):
    query, start_index, mapped_query = generate_query_sequence(nucleotide_array, probs_array)
    return query,start_index, mapped_query
    