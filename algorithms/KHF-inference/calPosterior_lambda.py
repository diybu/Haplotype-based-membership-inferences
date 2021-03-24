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
from scipy.special import comb

from sympy import Symbol, solve, symbols
from itertools import permutations





def calPosterior(dunique, Lenumerate, dtypeSymbol,dtype):
	ltype = Lenumerate[0].keys()
	ntype = len(ltype)
	dabsent = {} #store absent probability for each haplotype
	dpresent = {} #store in DB probability for each haplotype

	for dsolve in Lenumerate:
		for stype in ltype:
                        idx = int(str(eval('stype')).split('x')[-1])
                        pfreq = dtype[dtype.keys()[idx]]

			cnt = int(dsolve[stype])
			if int(cnt) == 0:
                                logp = 1000*math.log(pfreq)
				#logp = 0.0
				#for xx in ltype:
				#	if xx == stype:
				#		continue
				#	logp += int(dsolve[xx])*math.log(dtypeSymbol[str(eval('xx'))])
				try:
					dabsent[stype].append(logp)
				except:
					dabsent[stype] = [logp]
			else:
                                #nCk = math.factorial(1000) / math.factorial(cnt) /math.factorial(1000-cnt)
                                nCk = int(comb(1000, cnt, exact=True))
                                logp = nCk+cnt*math.log(pfreq)+(1000-cnt)*math.log(1-pfreq)
                                #logp = 0.0
                                #for xx in ltype:
                                # 	logp += int(dsolve[xx])*math.log(dtypeSymbol[str(eval('xx'))])
                                try:
                                        dpresent[stype].append(logp)
                                except:
                                        dpresent[stype] = [logp]


	#for each haplotype, calculate posterior probability
	dpos = {}
	for stype in ltype:
		print stype
		idx = int(str(eval('stype')).split('x')[-1])
		skey = dtype.keys()[idx]
		if stype not in dabsent:
			print stype, 'must in the database'
			dpos[skey] = 0.0
		elif stype not in dpresent:
			print stype, 'must be absent'
			dpos[skey] = 1e06
		else:
			#print dpresent
			#print dabsent
			#const = abs(int(np.median(dpresent[stype])))-500
                        pconsPre = 0
                        for xx in dpresent[stype]:
                            print 'dpresent xx', xx
                            try:
                                pconsPre += math.exp(xx)
                            except:
                                if xx < 1:
                                    pconsPre += -1e06
                                else:
                                    pconsPre += 1e06
                        pconsAbs = 0
                        for xx in dabsent[stype]:
                            print 'dabsent xx', xx
                            try:
                                pconsAbs += math.exp(xx)
                            except:
                                if xx < 1:
                                    pconsAbs += -1e06
                                else:
                                    pconsAbs += 1e06
			#pconsPre = sum([math.exp(xx) for xx in dpresent[stype]])
			#pconsAbs = sum([math.exp(xx) for xx in dabsent[stype]])
			isymbol = 'x'+str(idx)
			pPre = dtypeSymbol[str(eval('isymbol'))]
			pinDB = pconsPre*pPre/(pconsPre*pPre+pconsAbs*(1-pPre))
                        flambda = pconsPre*pPre/(pconsAbs*(1-pPre)+1e-06)
			dpos[skey] = flambda
			print pconsPre, pPre, pconsPre*pPre
			print pconsAbs, 1-pPre, pconsAbs*(1-pPre)
			print stype, 'the prob in DB', pinDB
			print stype, 'the freq in DB', pPre

	for stype in dunique.keys():
		print stype, 'must in the database'
		idx = int(str(eval('stype')).split('x')[-1])
                skey = dtype.keys()[idx]
		dpos[skey] = 0.0

	print 'dtype', dtype
	print 'dpos', dpos
	return dpos



if __name__ == '__main__':




	pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpCntSubset.pk', 'rb')
	dsnpFreqSub = pk.load(pkl_file)
	pkl_file.close()

	pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpCntWholeDB.pk', 'rb')
	dsnpFreq = pk.load(pkl_file)
	pkl_file.close()

	pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpAllele.pk', 'rb')
	dsnpAllele = pk.load(pkl_file)
	pkl_file.close()

	dcode = {'A':'1', 'C':'2', 'G':'3', 'T':'4'}
	ntotal = dsnpFreq['TOTAL']
	#ntotal = dsnpFreqSub['TOTAL']


	path = '/data2/diybu/beaconHaplopAttack/haploblock/sampling/solutionOuts_multiBlock/'

        for root, dirs, files in os.walk(path+'multiple/'):
                for fname in files:
			pkl_file = open(path+'multiple/'+fname, 'rb')
		        LL = pk.load(pkl_file)
		        pkl_file.close()
		
		        lblock = LL[0]
		        dunique = LL[1]
		        Lenumerate = LL[2]
		
		
		        lmarker = lblock[0] #idx of marker, start from 1
		        dtype = lblock[1] #type & freq A:1, C:2, G:3, T:4
		
		
			##obtain the prior probability
			dtypeSymbol = {}
			for i in range(len(dtype.keys())):
				symbol = 'x'+str(i)
				freq = dtype[dtype.keys()[i]]
				dtypeSymbol[symbol] = 1-math.exp(1000*math.log(1-freq))
		
			dpos = calPosterior(dunique, Lenumerate, dtypeSymbol, dtype)

			sout = fname.split('_')[-1]		
		        output = open('/data2/diybu/beaconHaplopAttack/haploblock/sampling/posteriorResult_lambda_multiBlock/'+sout, 'wb')
		        pk.dump([lblock,dpos], output)
		        output.close()


        for root, dirs, files in os.walk(path+'unique/'):
                for fname in files:                        
			pkl_file = open(path+'unique/'+fname, 'rb')
			LL = pk.load(pkl_file)
			pkl_file.close()

			lblock = LL[0]
			dunique = LL[1]
	
			dtype = lblock[1]	
			dpos = {}		
                        for i,v in dunique.items():
                                skey = dtype.keys()[int(str(i)[1:])]
                                if v>0:
                                    dpos[skey] = 0.0
                                else:
                                    print skey, 'not in DB'
                                    dpos[skey] = 1e-06



			sout = fname.split('_')[-1]
                        output = open('/data2/diybu/beaconHaplopAttack/haploblock/sampling/posteriorResult_lambda_multiBlock/'+sout, 'wb')
                        pk.dump([lblock,dpos], output)
                        output.close()
