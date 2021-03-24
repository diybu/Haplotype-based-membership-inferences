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
from scipy.stats import poisson

from sympy import Symbol, solve, symbols
from itertools import * 
from scipy.special import gamma
import operator
from more_itertools import locate


def multinomialLogPmf(lcnt, lprob):
	lpk = [lcnt[i]*math.log(lprob[i]) for i in range(len(lcnt))]
	fpmf = math.log(gamma(sum(lcnt+1))) - sum([math.log(gamma(xxx+1)) for xxx in lcnt]) + sum(lpk)
	
	return fpmf



def multinomialLogCdf(lcnt, lprob):
#P(X1<=n1,X2<=n2,...,Xk<=nk)
#approximation from "A Representation for Multinomial Cumulative Distribution Functions"
#Bruce Levin, The Annals of Statistics, v.9, n.5, pp.1123-1126, 1981

	#print 'lprob', lprob
	#print 'lcnt', lcnt	

	k = 0
	gamma1 = 0.0
	gamma2 = 0.0
	sum_s2 = 0.0
	sum_mu = 0.0
	sp = 0.0
	x = 0.0
	x2 = 0.0
	PWN = 0.0
	s = ntotal	
	log_cdf = -poisson.logpmf(ntotal, s)

	K = len(lcnt)
	for i in range(K):
		sp = s*lprob[i]
		pcdf = poisson.cdf(lcnt[i],  sp)
		log_cdf += poisson.logcdf(lcnt[i],  sp)
	
		mu = sp*(1 - poisson.pmf(lcnt[i], sp)/pcdf)
		s2 = mu - (lcnt[i] - mu)*(sp - mu)

		mr = lcnt[i]
		mf2 = sp * mu - mr * (sp - mu)

		mr *= lcnt[i] - 1
		mf3 = sp * mf2 - mr * (sp - mu)

		mr *= lcnt[i] - 2
		mf4 = sp * mf3 - mr * (sp - mu)

		m2 = mf2 + mu*(1-mu)
		mu3 = mf3 + mf2*(3-3*mu) + mu*(1+mu*(-3+2*mu))
		mu4 = mf4 + mf3*(6-4*mu) + mf2*(7+mu*(-12+6*mu))+mu*(1+mu*(-4+mu*(6-3*mu)))

		gamma1 += mu3
		gamma2 += mu4 - 3*s2*s2
		sum_mu += mu
		sum_s2 += s2

	sp = math.sqrt(sum_s2)
	gamma1 /= sum_s2 * sp
	gamma2 /= sum_s2 * sum_s2

	x = (ntotal-sum_mu)/sp
	x2 = x*x
	PWN = -x2/2 + math.log(1+gamma1/float(6)*x*(x2-3)+gamma2/float(24)*(x2*x2-6*x2+3)+gamma1*gamma1/72*(((x2-15)*x2+45)*x2-15)) - math.log(2*math.pi)/2  - math.log(sp)

	
	log_cdf += PWN
	

	return log_cdf


def postEquation(logpin, logpabs, logpPrior, logpPriorabs):

	pseudo = 10
        print 'logpin, logpabs, logpPrior, logpPriorabs', logpin, logpabs, logpPrior, logpPriorabs
        try:
                tmp1 = math.exp(logpin + logpPrior)
        except:
                if logpin + logpPrior < 0:
                        tmp1 = -1e06
                else:
                        tmp1 = 1e06
        try:
                tmp2 = math.exp(logpabs + logpPriorabs)
        except:
                if logpabs + logpPriorabs < 0:
                        tmp2 = -1e06
                else:
                        tmp2 = 1e06
        fpos = tmp1/(tmp1+ tmp2)
	return fpos


def mincnt(labsentIdx, ntype, nmustIn):
	lcnt = [ntotal-nmustIn for i in range(ntype)]
	for idx in labsentIdx:
		lcnt[idx] = 0

	return lcnt


def mincnt_in(labsentIdx, ntype, nmustIn):
        lcnt = [ntotal-1-nmustIn for i in range(ntype)]
        for idx in labsentIdx:
                lcnt[idx] = 0

        return lcnt


def mincnt_enum(llexist, ntype):
	nfree = ntotal - sum(llexist)
	lcnt = []
	for idx in range(ntype):
		if llexist[idx] == 0:
			lcnt.append(0)
		else:
			lcnt.append(nfree)	

	return lcnt



def checkComb(lcomb, lenumerate, lmustIn, labsentIdx):
#check if the combination meets the constraint 
	bpass = True
	for lidx in lenumerate:
		#print 'lidx', lidx
		#print 'lcomb', lcomb
		ltmp = [lcomb[i] for i in lidx]
		#print 'ltmp', ltmp		
		if 1 not in ltmp:
			bpass = False
			if bpass == False:
				#print 'lcomb without 1'
				break
		ltmp = [lcomb[i] for i in lmustIn]
		if 0 in ltmp:
			bpass = False
			if bpass == False:
				#print 'lmustIn has 0'
				break
		ltmp = [lcomb[i] for i in labsentIdx]
		if 1 in ltmp:
			bpass = False
			if bpass == False:
				#print 'labsentIdx has 1'
				break

	return bpass	



def calPosterior(Lsolve, ltypeSymbol):
        labsentIdx = Lsolve[0]
        lmustIn = Lsolve[1]
        lenumerate = Lsolve[2]

	if [] in lenumerate:
		lenumerate.remove([])
	if lenumerate == []:
		return []

	#print 'labsentIdx', labsentIdx
	#print 'lmustIn', lmustIn
	#print 'lenumerate', lenumerate

	print 'ltypeSymbol', ltypeSymbol #prior probability, frequency calculated from public DB


	ntype = len(ltypeSymbol)
	ltype = range(ntype)
	dabsent = {} #store absent probability for each haplotype
	dpresent = {} #store in DB probability for each haplotype

	lpos = [-1.0 for i in range(ntype)]
	for idx in labsentIdx:
		lpos[idx] = 0.0
	
	for idx in lmustIn:
		lpos[idx] = 1.0

	#print 'lpos', lpos

	#calculate the list of types w.o. constraint
	lsetEnu = []
	for l in lenumerate:
		for i in l:
			if i not in lsetEnu:
				lsetEnu.append(i)
	#print 'lsetEnu', lsetEnu
	lwo = list(set(ltype) - set(lmustIn) - set(labsentIdx) -set(lsetEnu))

	pcnt = 1e-06
	for idx in lwo:
		lcnt = mincnt(labsentIdx, ntype, len(lmustIn))
		lcnt[idx] = 0
		#print 'lcnt', lcnt
		logpab = multinomialLogCdf(lcnt, ltypeSymbol)
		#print 'logpabs', logpabs	

		lcnt_in = mincnt_in(labsentIdx, ntype, len(lmustIn))	
		logpin = multinomialLogCdf(lcnt_in, ltypeSymbol)
		#print 'logpin', logpin
		print 'ltypeSymbol[idx]', ltypeSymbol[idx], idx
                if ltypeSymbol[idx]==0.0:
                        lpos[idx] = 0.0
                elif ltypeSymbol[idx]==1.0:
                        lpos[idx] = 1.0
                else:
            		logpPriorin = math.log(ltypeSymbol[idx])
            		logpPriorabs = math.log(1 - ltypeSymbol[idx])
            		print 'logpin, logpabs, logpPriorin, logpPriorabs', logpin, logpabs, logpPriorin, logpPriorabs
            		fpos = postEquation(logpin, logpabs, logpPriorin, logpPriorabs)
            		#print 'fppost', fppost
            		lpos[idx] = fpos

	
	#enumerate the list of types among which one of them must be in DB
	denumerate = {}
	
	for l in lenumerate:
		for idx in l:
			if idx not in denumerate.keys():
				denumerate[idx] = set(l)
			else:
				denumerate[idx] = set(list(denumerate[idx]) + l)
				for j in set(denumerate[idx]-set([idx])):
					if j not in denumerate.keys():
						denumerate[j] = denumerate[idx] 
					else:
						denumerate[j] = set(list(denumerate[j]) + list(denumerate[idx]))
					
						
	#print 'denumerate, enumerate the list of types among which one of them must be in DB', denumerate



######################calculation of prob w/ constraint################################

	for idx in denumerate.keys():

		print 'working on idx, w/ constraint', idx	

		#idx not in DB, one of other snps w/ constraint must in DB: 1 - none of other snps w/ constraint in DB
		llexist0Idx = labsentIdx + list(denumerate[idx])
		lcnt_1 = mincnt(labsentIdx, ntype, len(lmustIn))
		lcnt_none_in = mincnt(llexist0Idx, ntype, len(lmustIn))
                lcnt_1[idx] = 0
		lcnt_none_in[idx] = 0
		
                #print 'lcnt', lcnt
                logpabs = multinomialLogCdf(lcnt_1, ltypeSymbol) - multinomialLogCdf(lcnt_none_in, ltypeSymbol)
                print 'logpabs', logpabs       

                lcnt_in = mincnt_in(labsentIdx, ntype, len(lmustIn))
                logpin = multinomialLogCdf(lcnt_in, ltypeSymbol)
                print 'logpin', logpin
                logpPriorin = math.log(ltypeSymbol[idx])
                logpPriorabs = math.log(1 - ltypeSymbol[idx])
                print 'logpin, logpabs, logpPriorin, logpPriorabs', logpin, logpabs, logpPriorin, logpPriorabs
                fpos = postEquation(logpin, logpabs, logpPriorin, logpPriorabs)
                print 'fpos', fpos
                lpos[idx] = fpos
	


			

	print 'lpos', lpos
	return lpos



if __name__ == '__main__':

        pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpCntWholeDB.pk', 'rb')
        dsnpFreq = pk.load(pkl_file)
        pkl_file.close()


	pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpCntSubset.pk', 'rb')
	dsnpFreqSub = pk.load(pkl_file)
	pkl_file.close()

	pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpAllele.pk', 'rb')
	dsnpAllele = pk.load(pkl_file)
	pkl_file.close()

	dcode = {'A':'1', 'C':'2', 'G':'3', 'T':'4'}
	ntotal = dsnpFreq['TOTAL']
	#ntotal = dsnpFreqSub['TOTAL']


	path = '/data2/diybu/beaconHaplopAttack/haploblock/beaconDB/solutionOuts_backup/'	

        for root, dirs, files in os.walk(path):
                for fname in files:
                        pkl_file = open(path+fname, 'rb')
                        LL = pk.load(pkl_file)
                        pkl_file.close()

                        lblock = LL[0]
			Lsolve = LL[1]

			#lmarker = lblock[0] #idx of marker, start from 1
			dtype = lblock[1]
			print 'dtype', dtype


                        ##obtain the prior probability
                        ltypeSymbol = []
                        for i in range(len(dtype.keys())):
                                symbol = i
                                freq = dtype[dtype.keys()[i]]
				print 'freq', freq
				fprior = 1-math.exp(1000*math.log(1-math.exp(freq)))
				#print 'fprior', fprior
                                ltypeSymbol.append(fprior)

			#print dtypeSymbol

			
			#dtype.values() -- frequency
			#ltypeSymbol = []
			#for isymbol in dtype.values():
			#	ltypeSymbol.append(2.71828**isymbol)
			
			lpos = calPosterior(Lsolve, ltypeSymbol)
			if lpos == []:
				continue

			dpos = {}
			for i in range(len(dtype.keys())):
				stype = dtype.keys()[i]
				dpos[stype] = lpos[i]

			sout = fname.split('_')[-1] 
			output = open('/data2/diybu/beaconHaplopAttack/haploblock/beaconDB/posteriorResult/beaconDB_'+sout, 'wb')
                        pk.dump([lblock,dpos], output)
                        output.close()


