#!/usr/bin/env python


# this script is to sample snp distribution based on haplotypes in one block
# by Diyue Bu


import sys
import os
import random
import time
import math
import pickle as pk
import numpy as np
from scipy.linalg import eigh, cholesky
from scipy.stats import norm

from sympy import Symbol, solve, symbols
from itertools import permutations



def pval(ltest, l_null):
        dThr = 0.95
        nthr = sorted(l_null)[int(len(l_null)*(1 - dThr))]
        pseudo = float(pow(10, -290))

        lpval = []
        for lrt in ltest:
                queryNorm = (lrt - mean(l_null))/(pseudo+std(l_null))
                pvalue = norm.cdf(queryNorm)
                lpval.append(pvalue)


        return lpval
		

def readFile(path):

        Llambda = []
        for root, dirs, files in os.walk(path):
                for fname in files:
                        pkl_file = open(path+fname, 'rb')
                        LL = pk.load(pkl_file)
                        pkl_file.close()

                        lblock = LL[0]
                        dpos = LL[1]
                        Llambda += dpos.values()
        return Llambda




if __name__ == '__main__':



	path_test = '/data2/diybu/beaconHaplopAttack/haploblock/beaconDB/posteriorResult/'
        path_null = '/data2/diybu/beaconHaplopAttack/haploblock/beaconDB/posteriorResult_null/'
        
        ltest = readFile(path_test)
        lnull = readFile(path_null)

        lpval = pval(ltest,lnull)

        print 'p-value', sum([True for x in lpval if x<=0.05 else False])

