#!/usr/bin/python
#bioinfo.py

__author__ = '''Hyunju Kim'''

import sys
import os
import random as ran
from math import log
from optparse import OptionParser, OptionGroup
from scipy import *
import matplotlib.pyplot as plt
import numpy as np
import itertools
from collections import defaultdict
import infodyn_entra
import operator
import draw_plots


def main(args):

    network_name = '%s-yeast-cell-cycle'%args[1]
    network_index = args[1]
    maxStep = int(args[2])
    historyLength = int(args[3])
    random_index = int(args[4])
    
    NODES_LIST_FILE = 'network-input/%s/%s-network-nodes-list.txt'%(network_name, network_name)
    nodes_list = []
    for line in open(NODES_LIST_FILE, 'r').readlines():
        items = [x.strip() for x in line.rstrip().split('\t')]
        if line[0] == '#':
            continue
        if line=='':
            continue
        nodes_list.append(items[0])


    input_file_name1 = 'time-series-random/%s-lr%d-step%d-trans0.dat'%(network_index, random_index, maxStep)
    input_file1 = open( input_file_name1, 'r')
    
    result_ai = open('results/ai-%s-lr%d-step%d-trans0-h%d.dat'%(network_index, random_index, maxStep, historyLength),'w')
    result_te = open('results/te-%s-lr%d-step%d-trans0-h%d.dat'%(network_index, random_index, maxStep, historyLength),'w')


    timeSeries = {}
    for n in nodes_list:
        timeSeries[n] = {}
        for si in range(0, pow(2,len(nodes_list))):
            timeSeries[n][si] = []

    for line in input_file1:
        if line[0] == '#':
            continue
        items = [x.strip() for x in line.rstrip().split('\t')]
        si = int(items[0])
        k = 3
        for n in nodes_list:
            timeSeries[n][si].append(int(items[k]))
            k +=1

    print 'AI'
    AI = {}
    for n in nodes_list:
        AI[n] = infodyn_entra.comAI(timeSeries[n], historyLength, nodes_list)
        result_ai.write('%s\t%f\n'%(n, AI[n]))
        print n, AI[n]
    print 'done AI'

    result_file_name = 'results/ai-scale-%s-lr%d-step%d-trans0-h%d.dat'%(network_index, random_index, maxStep, historyLength)
    viz_file_name = 'viz/ai-scale-%s-lr%d-step%d-trans0-h%d.pdf'%(network_index, random_index, maxStep, historyLength)
    draw_plots.rev_sorted_plot_file(AI, result_file_name, viz_file_name) ### plot and result file for AI scale



    print 'TE'
    TE =  defaultdict(float)
    for v in nodes_list:
        for n in nodes_list:
            TE[(v, n)] = infodyn_entra.comTE(timeSeries[v], timeSeries[n], historyLength, nodes_list)
            result_te.write('%s\t%s\t%f\n'%(v, n,TE[(v, n)] ))
            print v, n, TE[(v, n)]
    print 'done TE'


    TEdata = {}
    #for v in nodes_list:
    #for n in nodes_list:
    for v in nodes_list:
        for n in nodes_list:
            TEdata['%s-%s'%(v,n)] = TE[(v, n)]


    result_file_name = 'results/te-scale-%s-lr%d-step%d-trans0-h%d.dat'%(network_index, random_index, maxStep, historyLength)
    viz_file_name = 'viz/te-scale-%s-lr%d-step%d-trans0-h%d.pdf'%(network_index, random_index, maxStep, historyLength)
    draw_plots.rev_sorted_plot_file(TEdata, result_file_name, viz_file_name) ### plot and result file for TE scale



if __name__=='__main__':
    main(sys.argv)
