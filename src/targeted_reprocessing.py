#name: Mélanie Ghaby


################################### PART 3 ###################################
#Given a selected subarray, reprocesses it into a dictionary where for every
#word of a given length, any nucleotide with a probability less than 0.95 will 
#be tested out and so all possible words are recombining. the ones whose probability
#is less than 0.01 are rejected and not added to the dictionary
#uses concurrent programming to speed up the process

import itertools
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import json

#given a word and the respective probabilities of each nucleotide composing it
#it generates all possible words and their total probabilities by testing out
#every nucleotide option for the predicted nucleotides with less than a specific
#threshold of certainty
def generate_possible_words(examined_word, examined_probs, threshold, word_size):
    possible_words = []
    probabilities = []

    for i in range(word_size):
        nucleotide = examined_word[i]
        prob = examined_probs[i]
        
        #consider all nucleotides for a position whose
        if prob < threshold:
            possible_nucleotides = 'ACGT'  
            possible_words.append(possible_nucleotides)
            current_probabilities = []
            for nuc in ['A','C','G','T']:
                if nuc==nucleotide:
                    current_probabilities.append(prob)
                else:
                    current_probabilities.append((1-prob)/3)
                    
            probabilities.append(current_probabilities)  
        else: #keep nucleotide and probability as is if above threshold
            possible_words.append([nucleotide])  
            probabilities.append([prob])  
    
    #generating all possible combinations for nucleotides and their 
    #respective probabilities and converting them into list
    all_word_combinations = list(itertools.product(*possible_words))
    all_probabilities = list(itertools.product(*probabilities))

    #computing the probability of each word existing
    #every probability and combination are associated with zip
    word_probabilities = []
    for comb, prob_comb in zip(all_word_combinations, all_probabilities):
        word = ''.join(comb)
        prob= 1 #starting with 1 as a base
        for p in prob_comb:
            prob *= p  
        if prob > 0.01: #rejecting words whose probability falls below 0.01
            word_probabilities.append((word, prob))
        
    return word_probabilities

#saves the probabilisitc dictionary in a file
def save_dict_to_json(word_dict, word_size):
    filename = f'w_{word_size}_with_p.json'
    with open(filename, 'w') as json_file:
        json.dump(word_dict, json_file, indent=4)
    
    return filename

#building a probabilistic dictionary using parallel computing
def build_probabilistic_blast_database(nucleotide_array, prob_array, word_size, threshold):
    #initializing the dictionary
    word_dict = defaultdict(list)
    
    #analyszes a word starting at an index i of the subarray of nucleotides
    #and explores variations of nucleotides for the same word at positions
    #with probabilities that are too low
    def analyze_word(i):
        examined_word = nucleotide_array[i:i + word_size]
        examined_probs = prob_array[i:i + word_size]
        possible_words_with_probs = generate_possible_words(examined_word, examined_probs, threshold, word_size)
        
        #storing tuples of (generated word, index in database, probability)
        #using vectorizaiton to make things faster
        return [(word, i, prob) for word, prob in possible_words_with_probs]

    #processing all possible words in parallel with parallel computing
    with ThreadPoolExecutor() as executor:
        #analyze_word function executing in parallel many cores
        total_possible_words = executor.map(analyze_word, range(len(nucleotide_array) - word_size + 1))
    
    
    #adding the entries to the dictionary
    for result in total_possible_words:
        for word, index, prob in result:
            word_dict[word].append((index, prob))
    

    filename = save_dict_to_json(word_dict, word_size)
    
    return filename