# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 13:52:00 2018

@author: Jarnd
"""

import numpy as np
import math as mth
import itertools as itt
import copy

#X = np.array([[0, 1],[1, 0]]);
#Z = np.array([[1, 0],[ 0, -1]]);
#Y = 1j*np.dot(X,Z);
#I = np.matmul(X,X);

def bitstring_to_pauli(stab_list):
    # This function takes an numpy array of size 2*1 with every entry a n-bit 
    # string, representing a stabilizer, and returning a string of length n 
    # representing the stabilizer as its Pauli's
    pau_dict = {0:'I',1:'Z',2:'X',3:'Y'};
    stab_str_list = []
    if type(stab_list[0][0]) == np.int32:
        stab_syn = list(stab_list[0]+2*stab_list[1]);
        stab_str = ''
        for entry in stab_syn:
            stab_str += pau_dict[entry]
        return stab_str
    
    for stab in stab_list:
        stab_syn = list(stab[0]+2*stab[1]);
        stab_str = ''
        for entry in stab_syn:
            stab_str += pau_dict[entry]
        stab_str_list.append(stab_str)
    return stab_str_list

def comp_comm(A,B):
    # Returns the commutation relation between two stabilizers A and B represented as bitstrings.
    # Return 0 for commutation and 1 for anticommutation
    return (np.dot(A[0],B[1])+np.dot(A[1],B[0])) % 2

def stabgroup(gen_list):
    k = len(gen_list[0][0]);
    Z = np.zeros((k,), dtype=int)   # Defining an empty array seems not to work
    X = np.zeros((k,), dtype=int)
    for item in gen_list:
        Z = np.vstack((Z,item[0]))
        X = np.vstack((X,item[1]))
    stabgroup = (np.delete(Z,0,0), np.delete(X,0,0)) # Now I need to delete these, can I do this in another way?
    return stabgroup

def check_stabgroup(stab_group):
    # Checks is all stabilizers in the stabgroup commute
    ### TO DO: put in check if all provided elements are indeed generators (so no multiplications)
    ### Check rank of matrix
    Z = stab_group[0];
    X = stab_group[1];
    
    if np.count_nonzero((Z @ np.transpose(X) - X @ np.transpose(Z))%2) != 0:
        return [False, (Z @ np.transpose(X) - X @ np.transpose(Z))%2]
    elif rank(gen_group_to_matrix(stab_group)) != len(Z[:,0]):
        return [False, rank(gen_group_to_matrix(stab_group))]
    else:
        return True

def gen_matrix_to_stab_list(S):
    gen_list = []
    h,w = np.shape(S)
    for i in range(h):
        gen_list.append([S[i,:int(w/2)],S[i,int(w/2):]])
    return gen_list
        

def gen_group_to_matrix(gen_group):
    return np.concatenate((gen_group[0],gen_group[1]),axis=1)

def pauli_mult(P_A,P_B):
    # Computes the multiplication of two Pauli's as a bitstringsum.
    ## TO DO(?): Add phase?
    return (P_A+P_B) % 2

def list_stabgroup_elements(gen_group,stab, notation = None):
    element_list = []
    if notation is None:
        notation = 'B'
    for item in itt.product(range(2),repeat = len(gen_group[0])):
        element = copy.deepcopy(stab)
        element[0] = np.sum(gen_group[0]*np.array(item)[:, None],axis = 0) % 2;
        element[1] = np.sum(gen_group[1]*np.array(item)[:, None],axis = 0) % 2;
        if notation == 'B': element_list.append(element)
        elif notation == 'P': element_list.append(bitstring_to_pauli(element))
    return element_list

def swap_row(mat,row1,row2):
    mat[[row1,row2]] = mat[[row2,row1]];
    return mat

def ech(K):
    S = K.copy()
    rank = 0
    for i in range(np.shape(S)[1]):
        ones = np.where(S[rank:,i] == 1)[0]
        if len(ones) > 0:
            S = swap_row(S,ones[0]+rank,rank)
            ones[0] = rank
            for j in range(1,len(ones)):
                S[ones[j]+rank,:] = (S[ones[j]+rank,:] - S[ones[0],:])%2
            rank +=1;
    return S
             
def rank(K):
    S = ech(np.copy(K))
    r = np.shape(S)[0]
    for row in range(r):
        if not np.any(S[row,:]):
            r -= 1
    return r 

def piv_columns(K):
    S = ech(K.copy())
    piv_column = []
    for row in range(rank(S)):
        piv_column.append(np.where(S[row,:] == 1)[0][0])
    return piv_column

def red_ech(K):
    S = ech(K.copy())
    columns = piv_columns(S)
    for column in columns:
        ones_ind = list(np.where(S[:,column] == 1)[0])
        if len(ones_ind) > 1:
            for entry in range(0,len(ones_ind)-1):
                S[ones_ind[entry],:] = (S[ones_ind[entry],:]+S[ones_ind[-1],:])%2
    return S

def swap_XZ(vect):
    l = np.size(vect)
    h = int(l/2)
    vectn = np.zeros(l)
    vectn[:h] = vect[h:]
    vectn[h:] = vect[:h]
    return vectn.astype(int)

def nullspace_basis(K):
    S = red_ech(K.copy()) 
    w,l = np.shape(K)
    basis = []
    non_piv = [x for x in [x for x in range(l)] if x not in piv_columns(K)]
    for column in non_piv:
        vect = np.zeros(l)
        vect[:w] = S[:w,column]
        vect[column] = 1
        basis.append(swap_XZ(vect.astype(int)))
        
    return basis

def logical_pauli(K):
    S = red_ech(K.copy())
    paulis = []
    basis = nullspace_basis(S)
    for vect in basis:
        newmatr = np.vstack((K,vect))
        if rank(newmatr) != rank(K):
            paulis.append(vect)
    return paulis
    
        
#def list_detectable_errors(gen_group):
#    S = np.concatenate((gen_group[0],gen_group[1]),axis=1);
        
                
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        