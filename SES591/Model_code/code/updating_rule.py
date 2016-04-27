#!/usr/bin/python
#attractor-analysis.py
#last update: 17 DEC 2015

__author__ = '''Hyunju Kim'''


import os
import sys
import numpy as np
import networkx as nx
from collections import OrderedDict

import input_net as inet



def sigmoid_updating(net, prevState):

	currState = {}
	#### compute the current states of nodes in the net ####

	for u,v in prevState.iteritems():
		#### CycD
		if prevState['CycD'] == 1:
			currState['CycD'] = 1
		else:
			currState['CycD'] = 0
		#### Rb
		if (prevState['CycD'] == 0 and prevState['CycE'] == 0 and prevState['CycA'] == 0 and prevState['CycB'] == 0) or (prevState['P27'] == 1 and prevState['CycD'] == 0 and prevState['CycB'] == 0):
			currState['Rb'] = 1
		else:
			currState['Rb'] = 0
		#### E2F
		if (prevState['Rb'] == 0 and prevState['CycA'] == 0 and prevState['CycB'] == 0) or (prevState['P27'] == 1 and prevState['Rb'] == 0 and prevState['CycB'] == 0):
			currState['E2F'] = 1
		else:
			currState['E2F'] = 0
		#### CycE
		if prevState['E2F'] == 1 and prevState['Rb'] == 0:
			currState['CycE']= 1
		else: 
			currState['CycE'] = 0
		#### CycA
		if (prevState['E2F'] == 1 and prevState['Rb'] == 0 and prevState['Cdc20'] == 0 and (prevState['Cdh1'] == 0 or prevState['UbcH10'] == 0)) or (prevState['CycA'] == 1 and prevState['Rb'] == 0 and prevState['Cdc20'] == 0 and (prevState['Cdh1'] == 0 or prevState['UbcH10'] == 0)):
			currState['CycA'] = 1
		else:
			currState['CycA'] = 0
		#### P27
        if (prevState['CycD'] == 0 and prevState['CycE'] == 0 and prevState['CycA'] == 0 and prevState['CycB'] == 0) or (prevState['P27'] == 1 and (prevState['CycE'] == 0 or prevState['CycA'] == 0) and prevState['CycB'] == 0 and prevState['CycD'] == 0):
            currState['P27'] = 1
        else:
			currState['P27'] = 0
		#### Cdc20
        if prevState['CycB'] == 1:
			currState['Cdc20'] = 1
        else:
			currState['Cdc20'] = 0
		#### Cdh1
        if (prevState['CycA'] == 0 and prevState['CycB'] == 0) or (prevState['Cdc20'] == 1) or (prevState['P27'] ==1  and prevState['CycB'] == 0):
			currState['Cdh1'] = 1
        else:
			currState['Cdh1'] = 0
		#### UbcH10
        if prevState['Cdh1'] == 0 or (prevState['Cdh1'] == 1 and prevState['UbcH10'] == 1 and (prevState['Cdc20'] == 1 or prevState['CycA'] == 1 or prevState['CycB'] == 1)):
            currState['UbcH10'] = 1
        else:
			currState['UbcH10'] = 0
		#### CycB
        if prevState['Cdc20'] == 0 and prevState['Cdh1'] == 0:
			currState['CycB'] = 1
        else:
			currState['CycB'] = 0
			
	return currState
		
		################# begin: sigmoid_updating ######################

    
################# end: sigmoid_updating ########################

def main():
	print "updating_rule module is the main code."
	EDGE_FILE = '../data/example/example-net-edges.dat'
	NODE_FILE = '../data/example/example-net-nodes.dat'
	
	net = inet.read_network_from_file(EDGE_FILE, NODE_FILE)
	nodes_list = inet.build_nodes_list(NODE_FILE)
	
	#prevState = {'CycD':0.0, 'b':0.0, 'c':1.0}
	prevState = {}
	prevState['CycD'] = float(sys.argv[1])
	prevState['Rb'] = float(sys.argv[2])
	prevState['E2F'] = float(sys.argv[3])
	prevState['CycE'] = float(sys.argv[4])
	prevState['CycA'] = float(sys.argv[5])
	prevState['P27'] = float(sys.argv[6])
	prevState['Cdc20'] = float(sys.argv[7])
	prevState['Cdh1'] = float(sys.argv[8])
	prevState['UbcH10'] = float(sys.argv[9])
	prevState['CycB'] = float(sys.argv[10])
	#print "network state @ previous step", OrderedDict(sorted(prevState.items(), key=lambda t: t[0]))
	for v in nodes_list:
		print prevState[v],
	currState = sigmoid_updating(net, prevState)
	print '\n'
	#print "network state @ current step", OrderedDict(sorted(currState.items(), key=lambda t: t[0]))
	for v in nodes_list:
		print currState[v],

if __name__=='__main__':
    main()
