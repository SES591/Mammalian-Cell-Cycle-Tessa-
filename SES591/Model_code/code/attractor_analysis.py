#!/usr/bin/python
#bioinfo.py

__author__ = '''Hyunju Kim'''

import os
import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import input_net as inet
import updating_rule as ur
import math
import zlib

def binary_to_decimal(nodes_list, biStates, Nbr_States=2):  # more left in the nodes list means higher order of 2 in binary
    decState = 0
    x = len(nodes_list) -1
    for u in nodes_list:
        decState = decState + biStates[u]  * np.power(Nbr_States, x)
        x = x - 1
    return decState
    
def decimal_to_binary(nodes_list, decState, Nbr_States=2): # more left in the nodes list means higher order of 2 in binary
    biStates = {}
    x = len(nodes_list) -1
    for u in nodes_list:
        biStates[u] = decState / np.power(Nbr_States, x)
        decState = decState % np.power(Nbr_States, x)
        x = x - 1
        
    return biStates
def find_attractor(decStateTransMap):
    
    '''
        Arguments:
        1. decStateTransMap
        Return:
        1. attractor
    '''
    attractor_list = nx.simple_cycles(decStateTransMap) #in case of deterministic system, any cycle without considering edge direction will be directed cycle.
    attractors = {}
    attractors['fixed'] = []
    attractors['cycle'] = []

    for u in attractor_list:
        if len(u) == 1:
            attractors['fixed'].append(u)
        else:
            attractors['cycle'].append(u)

    return attractors

def net_state_transition(net, nodes_list, Nbr_States=2):

    '''
    Arguments:
               1. net
               2. Nbr_States
    Return:
               1. decStateTransMap
    '''
    
    Nbr_Nodes = len(net.nodes())
    Nbr_All_Initial_States = np.power(Nbr_States, Nbr_Nodes)
    
    decStateTransMap = nx.DiGraph()
    for prevDecState in range(0, Nbr_All_Initial_States):
        prevBiState = decimal_to_binary(nodes_list, prevDecState, Nbr_States)
        currBiState = ur.sigmoid_updating(net, prevBiState)
        currDecState = binary_to_decimal(nodes_list, currBiState, Nbr_States)
        decStateTransMap.add_edge(prevDecState, currDecState)
    return decStateTransMap
'''
def	biological_sequence(net, nodes_list, bio_initStates, fileName, Nbr_States=2):
    bioSeq = []
    currBiStates = bio_initStates
    finished = False
    while(not finished):
        oneDiff = 0
        prevBiStates = currBiStates.copy()
        bioSeq.append(prevBiStates)
        currBiStates = ur.sigmoid_updating(net, prevBiStates)
        for u in nodes_list:
            if abs(prevBiStates[u] - currBiStates[u]) > 0:
                oneDiff = 1
                break
        finished = (oneDiff < 1)

    OUTPUT_FILE  = open(fileName, 'w')
    OUTPUT_FILE.write('time step')
    for u in nodes_list:
        OUTPUT_FILE.write('\t%s'%(u))
    OUTPUT_FILE.write('\n')

    for i in range(len(bioSeq)):
        OUTPUT_FILE.write('%d'%i)
        for u in nodes_list:
            OUTPUT_FILE.write('\t%d'%(bioSeq[i][u]))
        OUTPUT_FILE.write('\n')
        return BioSeq
'''
def main():
	
	EDGE_FILE = '../data/fission-net/fission-net-edges.txt'
	NODE_FILE = '../data/fission-net/fission-net-nodes.txt'
	BIO_INIT_FILE = '../data/fission-net/fission-net-bioSeq-initial.txt'
	
	net = inet.read_network_from_file(EDGE_FILE, NODE_FILE)
	nodes_list = inet.build_nodes_list(NODE_FILE)
	bio_initStates = inet.read_init_from_file(BIO_INIT_FILE)
	decStateTransMap = net_state_transition(net, nodes_list)
	attractors = find_attractor(decStateTransMap)
	print attractors



    

if __name__=='__main__':
    main()
    
    
