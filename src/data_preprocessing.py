#name: Mélanie Ghaby

############################# PART 1 ##########################################
#This is where the initial processing happens.
#The probabilistic nucleotide array is put into a dictionary
#Of Words of some size K between [6,11]

import json

#importing the sequence of nucleotides
def nucleotides_to_array(filename):
    with open(filename, 'r') as file:
        sequence = file.read().replace('\n', '')   
    nucleotide_array = list(sequence)
    return nucleotide_array


#importing the probabilities of nucleotides of the sequence
def probs_to_array(filename):
    with open(filename, 'r') as file:
        numbers = file.read().split()
    number_array = [float(num) for num in numbers]
    return number_array

  
#making a dictionary out of nucleotides_array for a given word_size
def preprocess_nucleotide_array(nucleotide_array, word_size, output_filename):
    words_dictionary = {}

    #extracting the words of a given size
    for i in range(len(nucleotide_array) - word_size + 1):  
        word = nucleotide_array[i:i + word_size]  
        word_str = ''.join(word) 
        #adding word if not in dictionary already
        if word_str not in words_dictionary:
            words_dictionary[word_str] = []
        #adding corresponding position in the database as a value in dictionmary
        words_dictionary[word_str].append(i)

    with open(output_filename, 'w') as json_file:
        json.dump(words_dictionary, json_file, indent=4)

def data_preprocessor(word_size):
    #importing the nucleotide sequence and the sequence of probabilities
    #that were downloaded
    filename = 'predicted_database.txt'
    nucleotide_array = nucleotides_to_array(filename)
    filename = 'predicted_database_probs.txt'
    probs_array = probs_to_array(filename)

    #ensuring both length matches
    if len(nucleotide_array) != len(probs_array):
        print("error in data: lengths do not match")
      
    #making the dictionary according to the chosen word_size
    #preprocess_nucleotide_array(nucleotide_array, word_size, f"w_{word_size}.json")
    #information_about_sequences(nucleotide_array, probs_array)
    
    return nucleotide_array, probs_array


#################used for the report part#############
def information_about_sequences(nucleotide_array, probs_array):
    min_value = min(probs_array)
    max_value = max(probs_array)

    below_threshold_count = 0
    for p in probs_array:
        if p < 0.95:
            below_threshold_count += 1

    percentage_below_threshold = (below_threshold_count / len(probs_array)) * 100
    
    certitude_count = 0
    for p in probs_array:
        if p == 1.0:
            certitude_count += 1

    percentage_below_threshold = (below_threshold_count / len(probs_array)) * 100
    percenteage_of_certitude = (certitude_count/len(probs_array)) *100

    print(f"Least likely predicted nuc: {min_value}")
    print(f"Percentage of certain nucleotides: {percenteage_of_certitude:.2f}%")
    print(f"Percentage of values below 0.95: {percentage_below_threshold:.2f}%")

