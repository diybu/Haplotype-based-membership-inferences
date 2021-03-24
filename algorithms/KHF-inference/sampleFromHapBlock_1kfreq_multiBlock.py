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






def buildEquation(lsnpFreq, lalt, dtype):

	ltypes = dtype.keys()
	ntype = len(ltypes)
	#print 'number of haplotypes', ntype
	nsnp = len(lsnpFreq)
	#print 'number of snps in the block, length of vector b', nsnp

	tsymbols = symbols('x:'+str(ntype), positive=True)

	lequation = []
	lexist = []
	for j in range(nsnp):
		etmp = 0
		for i in range(ntype):
			stype = ltypes[i]
			if stype[j] == lalt[j]:
				etmp += tsymbols[i]
		if etmp not in lexist:
			lexist.append(etmp)
			lequation.append(etmp)
			#lequation[-1] -= round(lsnpFreq[j],2)
			#lequation[-1] -=lsnpFreq[j]
			lequation[-1] -= lsnpFreq[j]*ntotal

	#lequation.append(-1)
	#lequation.append(-sum(dtype.values()))
	lequation.append(-sum(dtype.values())*ntotal)
	for i in range(ntype):
		lequation[-1] += tsymbols[i] 				
	try:
		dsolve = solve(lequation, list(tsymbols), quick=True)
	except:
		dsolve = []

	bsolve = False
	if dsolve == []:
		print 'no solution'
                print lalt
                print lequation
                print lsnpFreq
                for _ in ltypes:
                        print _, dtype[_]*ntotal
		pass


	else:
		#bsolve = True
	        #print lalt
	        #print lequation
	        #print lsnpFreq
		#print dsolve
	        #for _ in ltypes:
	        #        print _, dtype[_]*ntotal
		if len(dsolve.keys()) != ntype:		
			bsolve = True
			print set(tsymbols) - set(dsolve.keys())	
			print dsolve
	                #print lalt
	                #print lequation
	                #print lsnpFreq
	                #for _ in ltypes:
	                #        print _, dtype[_]
			pass
		else:
			#print dsolve
			pass


	return dsolve, bsolve





def multiSolution(dsolve, nvar):
	dunique = {}
	Lenumerate = []

	print 'multiple solutions, enumerate all the results'
	tsymbols = symbols('x:'+str(nvar))
	lvar =  list(set(tsymbols) - set(dsolve.keys())) #store dependent var

	if len(lvar) >= 3:
		return {},[]

	#print 'number of vars', nvar
	#print 'independent vars', tvar
	for svar in dsolve.keys():
		solution = dsolve[svar]
                print 'svar, solution', svar, solution
		solution = int(round(solution,0))
		dunique[svar] = solution  #vars w/ unique solution
	ldependent = list(set(dsolve.keys()) - set(dunique.keys()))



	lpermutate = list(permutations(range(ntotal+1), r=len(lvar)))



	LEnumSolve = []
	nidx = len(ldependent)
	niter = len(lpermutate)
	for i in range(niter):
		perT = lpermutate[0]
		lpermutate.remove(perT)

		lequa = []
                for sdep in ldependent:
                        lequa.append(dsolve[sdep] - sdep)

                for j in range(len(lvar)):
                        lequa.append(lvar[j]-perT[j])

		print perT	
		print 'NEW PERMUTATION RESULT', lequa
		dsolveTemp = solve(lequa, lvar+ldependent, force=True)		
		print dsolveTemp
		if dsolveTemp == []:
			continue
		bcheck = all(j>=0 and j<=ntotal for j in dsolveTemp.values())
		if bcheck:
			LEnumSolve.append(dsolveTemp)
			print dsolveTemp
		else:
			pass

	del lpermutate

	for nrep in range(ntotal+1):

                lequa = []
                for sdep in ldependent:
                        lequa.append(dsolve[sdep] - sdep)

                for j in range(len(lvar)):
                        lequa.append(lvar[j]-nrep)

		print 'repeat permutation'
                print 'NEW PERMUTATION RESULT', lequa
                dsolveTemp = solve(lequa, lvar+ldependent, force=True) 
                print dsolveTemp
                if dsolveTemp == []:
                        continue
                bcheck = all(j>=0 and j<=ntotal for j in dsolveTemp.values())
                if bcheck:
                        LEnumSolve.append(dsolveTemp)
                        print dsolveTemp
                else:
                        pass

	
	return dunique, LEnumSolve
	
	#print '# of permutations', len(lpermutate)



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
	#ntotal = dsnpFreq['TOTAL']
	ntotal = dsnpFreqSub['TOTAL']

	del dsnpFreqSub

	cnt = 0
	lsolve = []
        for filename in os.listdir('/data2/diybu/beaconHaplopAttack/haploblock/hapblockTarget_multiBlock'):

                pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/hapblockTarget_multiBlock/'+filename, 'rb')
                lblock = pk.load(pkl_file)
                pkl_file.close()

		cnt += 1

                lmarker = lblock[0] #idx of marker, start from 1
                dprior = lblock[1] #type & prior A:1, C:2, G:3, T:4
                dfreq = lblock[2] #type & freq A:1, C:2, G:3, T:4
                fname = lblock[-1].split('.')[0]


		
	
		#get rsid
		lrsidx = []
		with open('/data2/diybu/beaconHaplopAttack/haploblock/genPed/'+fname+'.info', 'r') as f:
			for l in f:
				sid = l.strip()
				lrsidx.append(sid)


		lskey = []
		lsnpFreq = [0.0 for i in range(len(lmarker))]
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
			sref = dcode[lallele[0]]
			#lsnpFreq.append(dsnpFreq[skey]) #store the minor allele freq for each snp
			#print skey, dsnpAllele[skey], dsnpFreq[skey], dsnpFreqSub[skey]
	
			#sotre the minor allele freq for each snp for iterating each haplotype
			for stype in dfreq.keys():
				#print salt	
				if stype[i] == salt:	
					lsnpFreq[i] += dfreq[stype] 	
				elif stype[i] == sref:
					pass
				else:
					print 'stype, sref, salt do not match', idx, sid, sref, salt, stype[i]			

		#print '\t'.join(lskey)
                
                print 'begin equation'
		dsolve, bsolve = buildEquation(lsnpFreq, lalt, dfreq)
		if dsolve == []:
			print 'no solution for this block', fname
			continue
		print 'No of block', cnt
		print 'if multiple solutions', bsolve




		if bsolve:
                        print 'begin multisolution'
			dunique, Lenumerate = multiSolution(dsolve, len(dfreq.keys()))
                        print 'end multisolution'
			if dunique == {} and Lenumerate == []:
				continue

		        output = open('solutionOuts_multiBlock/multiple/LmultipleSolutions_block'+str(cnt)+'.pk', 'wb')
		        pk.dump([lblock, dunique, Lenumerate], output)
		        output.close()
			
			del Lenumerate	
			del dunique
			del dsolve
		
		else:

                        output_sub = open('solutionOuts_multiBlock/unique/LuniqueSolution_block'+str(cnt)+'.pk', 'wb')
                        pk.dump([lblock, dsolve], output_sub)
                        output_sub.close()
			del dsolve			


			#print 'all haoplotypes must be in the DB'
			#for stype in dtype.keys():
			#	print stype	
		lsolve.append(bsolve)



	print 'number of blocks having multiple solution', sum(lsolve)
	print 'number of blocks having unique solution', len(lsolve)-sum(lsolve)

	sys.exit()

	

