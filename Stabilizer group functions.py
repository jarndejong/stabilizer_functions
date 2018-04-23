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

def bitstring_to_pauli(stab):
    # This function takes an numpy array of size 2*1 with every entry a n-bit 
    # string, representing a stabilizer, and returning a string of length n 
    # representing the stabilizer as its Pauli's
    stab_syn = list(stab[0]+2*stab[1]);
    stab_str = ''
    pau_dict = {0:'I',1:'Z',2:'X',3:'Y'};
    for entry in stab_syn:
        stab_str += pau_dict[entry]
    return stab_str

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
        return False
    else:
        return True

def pauli_mult(P_A,P_B):
    # Computes the multiplication of two Pauli's as a bitstringsum.
    ## TO DO(?): Add phase?
    return (P_A+P_B) % 2

def list_stabgroup_elements(gen_group,stab):
    element_list = []
    for item in itt.product(range(2),repeat = 4):
        element = copy.deepcopy(stab)
        element[0] = np.sum(gen_group[0]*np.array(item)[:, None],axis = 0) % 2;
        element[1] = np.sum(gen_group[1]*np.array(item)[:, None],axis = 0) % 2;
        element_list.append(bitstring_to_pauli(element))
    return element_list

def list_detectable_errors(gen_group):
    S = np.concatenate((gen_group[0],gen_group[1]),axis=1);
    
       

stab1 = np.array([[0,1,1,0,0],[1,0,0,1,0]]);
stab2 = np.array([[0,0,1,1,0],[0,1,0,0,1]]);
stab3 = np.array([[0,0,0,1,1],[1,0,1,0,0]]);
stab4 = np.array([[1,0,0,0,1],[0,1,0,1,0]]);

gen_list = (stab1,stab2,stab3,stab4);
gen_group = stabgroup(gen_list)
print(check_stabgroup(gen_group))
print(list_stabgroup_elements(gen_group,stab1))