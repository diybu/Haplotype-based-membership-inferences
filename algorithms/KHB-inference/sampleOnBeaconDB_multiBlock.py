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






def constraint_haptype(s):
#remove the target haptypes in the sample to obtain the null distribution (haptype not in DB)

	l = list(filter(lambda x:x!=s, lhaptype))

	return l


def buildEquation(lsnpReturn, lalt, dtype):

	ltypes = dtype.keys()
	ntype = len(ltypes)
	print 'number of haplotypes', ntype
	nsnp = len(lsnpReturn)
	print 'number of snps in the block, length of vector b', nsnp

	tsymbols = symbols('x:'+str(ntype))

	lbysnp = [[] for i in range(nsnp)]
	for j in range(nsnp):
		for i in range(ntype):
			stype = ltypes[i]
			if stype[j] == lalt[j]:
				lbysnp[j].append(True)
			else:
				lbysnp[j].append(False)

	#print lsnpReturn
	labsentIdx = [] #store the index must not be in DB
	lpresentIdx = [] #store the index of types may be in DB	
	lfalse = [i for i,x in enumerate(lsnpReturn) if x==False]
	for idx in lfalse:
		ltfTypes = lbysnp[idx]
		ltmp = [i for i,x in enumerate(ltfTypes) if x==True] #types w/ the snp that is not in the DB must not present in DB	
		#print 'ltfTypes', ltfTypes
			
		labsentIdx += ltmp

	labsentIdx = list(set(labsentIdx))
	print 'labsentIdx', labsentIdx
		
	ltrue = sorted(list(set(range(len(lsnpReturn)))-set(lfalse)))
	print 'ltrue', ltrue
	lmustIn = []
	dmustIn = {i:False for i in range(ntype)}
	for idx in ltrue:
		ltfTypes = lbysnp[idx]
		#print 'ltfTypes', ltfTypes
		ltmp = [i for i,x in enumerate(ltfTypes) if x==True] #one of the types having this snp has to be present in DB
		if len(ltmp) == 1:
			if ltmp[0] not in lmustIn:
				lmustIn.append(ltmp[0])
				dmustIn[ltmp[0]] = True
			continue
		if ltmp not in lpresentIdx:
			lpresentIdx.append(ltmp)
		else:
			print ltmp
	print 'lmustIn', lmustIn
	print 'dmustIn', dmustIn
	print 'lpresentIdx', lpresentIdx

	lenumerate = []
	for ll in lpresentIdx:
		benumerate = True
		for idx in ll:
			if dmustIn[idx]:
				benumerate = False
				break
		if benumerate:
			lenumerate.append(ll)

	print 'lenumerate', lenumerate



	return [labsentIdx, lmustIn, lenumerate]






if __name__ == '__main__':



        pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpCntSubset.pk', 'rb')
        dsnpFreqSub = pk.load(pkl_file)
        pkl_file.close()



	pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpAllele.pk', 'rb')
	dsnpAllele = pk.load(pkl_file)
	pkl_file.close()

	dcode = {'A':'1', 'C':'2', 'G':'3', 'T':'4'}
	#ntotal = dsnpFreq['TOTAL']
	#ntotal = dsnpFreqSub['TOTAL']


	cnt = 0
	lsolve = []

        for filename in os.listdir('/data2/diybu/beaconHaplopAttack/haploblock/hapblockTarget_multiBlock'):

                pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/hapblockTarget_multiBlock/'+filename, 'rb')
                lblock = pk.load(pkl_file)
                pkl_file.close()

		cnt += 1

                lmarker = lblock[0] #idx of marker, start from 1
                dtype = lblock[1] #type & freq A:1, C:2, G:3, T:4
                fname = lblock[-1].split('.')[0]

		
	
		#get rsid
		lrsidx = []
		with open('/data2/diybu/beaconHaplopAttack/haploblock/genPed/'+fname+'.info', 'r') as f:
			for l in f:
				sid = l.strip()
				lrsidx.append(sid)


		lskey = []
		lsnpReturn = []
		lalt = []
		for i in range(len(lmarker)):
			idx = lmarker[i]
			idx = idx - 1
			sid = lrsidx[idx]
			lid = ['10'] + sid.split()[::-1]
			skey = ' '.join(lid)
			lskey.append(skey) #store the key to get the corresponding freq for each snp
			lallele = dsnpAllele[skey]
			salt = dcode[lallele[-1]]
			lalt.append(salt)
                        if skey in dsnpFreqSub:
                                lsnpReturn.append(True)
                        else:
                                lsnpReturn.append(False)  #store the minor allele existence for each snp
			#print skey, dsnpAllele[skey], dbReturn[skey]
	
		print 'length of the haploblock', len(lmarker)


		#print '\t'.join(lskey)	
		Lsolve = buildEquation(lsnpReturn, lalt, dtype)
	
                output = open('solutionOuts_multiblock/LbeaconSolutions_block'+str(cnt)+'.pk', 'wb')
                pk.dump([lblock, Lsolve], output)
                output.close()		

		#sys.exit()		





	

