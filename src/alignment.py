#name: Mélanie Ghaby

import numpy as np

################################### PART 5 ###################################
#Aligns the query sequence with the predicted range of belonging in the database
#Same formatting as my assignment 1

#gap penalty
def gap_penalty(prev_nuc_x):
    current_score = prev_nuc_x+np.log(0.002)
    return current_score
   
#initializes first row and first column
def initialization_x(query,dtb):
    matrix_x = np.zeros((len(query)+1, len(dtb)+1))
    matrix_x[0,0]=0
    
    #first column intialization
    for i in range(1, len(query)+1):
        matrix_x[i,0] = gap_penalty(matrix_x[i-1,0])
    #first row initializaton            
    for j in range(1, len(dtb)+1):
        matrix_x[0,j] = 0

    return matrix_x

#comptes score if adding both nucleotides to the sequence
#takes as input the current DP matrix filled out, current positions
#as well as the integers necessary for the soring
#outputs the score according to adding one nucleotide from each sequence
#mismatch whe probs subarray is 1 will be 0.002 bc < 0.003 = 0.01/3
def gapless_choice(matrix_x, i, j, query, dtb, probs_subarray):

    if query[i-1]==dtb[j-1]:
        return (matrix_x[i-1][j-1] + np.log(probs_subarray[j-1]))
    else:
        unmatched_probs = (1-probs_subarray[j-1])/3
        score_to_add = max(unmatched_probs, 0.002)
        return (matrix_x[i-1][j-1] + np.log(score_to_add))
 
#computes score if only adding a nucleotide to the sequence
def row_gap(matrix_x,i,j):
    return gap_penalty(matrix_x[i,j-1])

#computes score if only adding a nucleotide to the sequence
def col_gap(matrix_x, i,j):
    return gap_penalty(matrix_x[i-1,j])



#fills out one entry of the DP matrix in order to maximize it
#returns the max score, as well as a matrix called pointers_x
#pointers_x is a matrix storing the pointers
def matrix_maximizer(matrix_x,i,j,query,dtb,pointers_x,probs_subarray):
    diago = gapless_choice(matrix_x, i, j, query, dtb, probs_subarray)
    gapped_row = row_gap(matrix_x,i,j)
    gapped_col = col_gap(matrix_x, i,j)
    possibilities = [diago, gapped_row, gapped_col]
    max_score= max(possibilities)

    #storing the pointers
    if max_score == diago:
        pointers_x[i,j]=1
    elif max_score == gapped_row:
        pointers_x[i,j]=2
    else:
        pointers_x[i,j] =3
        
    return max_score, pointers_x


#fills out the DP matrix according to the matrix maximizer functoin
def DP_matrix_filling(query,dtb,probs_subarray):
    matrix_x = initialization_x(query,dtb)
    pointers_x = np.zeros((len(query)+1, len(dtb)+1))
    pointers_x[0,0] = 0
    for i in range(1, len(query)+1):
        pointers_x[i,0]= 3
    for j in range (1, len(dtb)+1):
        pointers_x[0,j] = 2
    last_query_gap=-1
    last_dtb_gap=-1
    
    for i in range(1, len(query) + 1):
        for j in range(1, len(dtb) + 1):
            maximizer_result = matrix_maximizer(matrix_x, i,j,query,dtb,pointers_x, probs_subarray)
            matrix_x[i, j] = maximizer_result[0]
            pointers_x = maximizer_result[1]
    
            
    return matrix_x, pointers_x



#backtracking from the pointers matrix in order to recover the alignment
def backtrack_alignment(query, dtb, pointers_x, matrix_x):
    aligned_query = []
    aligned_dtb = []
    i = len(query)
    #backtracking for the max in the last row
    j = np.argmax(matrix_x[-1])

    while i > 0 or j > 0:
        if pointers_x[i, j] == 1:  
            aligned_query.append(query[i - 1])
            aligned_dtb.append(dtb[j - 1])
            i -= 1
            j -= 1
        elif pointers_x[i, j] == 2:  
            aligned_query.append('-')
            aligned_dtb.append(dtb[j - 1])
            j -= 1
        elif pointers_x[i, j] == 3:  
            aligned_query.append(query[i - 1])
            aligned_dtb.append('-')
            i -= 1

    aligned_query.reverse()
    aligned_dtb.reverse()

    return ''.join(aligned_query), ''.join(aligned_dtb)


def full_algorithm(query_seq,database_seq,probs_subarray):
    query=query_seq
    dtb=database_seq
    DP_matrix_and_pointers = DP_matrix_filling(query,dtb,probs_subarray)
    DP_matrix = DP_matrix_and_pointers[0]
    pointers_matrix = DP_matrix_and_pointers[1]
    col_opt_score =  np.argmax(DP_matrix[-1])
    aligned_query, aligned_dtb = backtrack_alignment(query, dtb, pointers_matrix, DP_matrix)

    return col_opt_score,aligned_query,aligned_dtb
