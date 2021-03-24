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

from itertools import permutations
from pulp import *





def constraint_haptype(s):
#remove the target haptypes in the sample to obtain the null distribution (haptype not in DB)

	l = list(filter(lambda x:x!=s, lhaptype))

	return l



    
def findNewType(lknownType, lsnpFreq, lalt, dtype, nmiss):
##suppose there is one new haplotype in the DB, use integer linear programming to find this new type


	ltypes = dtype.keys()
	ntype = len(ltypes)
	print 'number of haplotypes', ntype+2
	for i in range(ntype):
		print 'real number of type', i, int(round(dtype.values()[i]*ntotal, 0))
	nsnp = len(lsnpFreq)
	print 'number of snps in the block, length of vector b', nsnp

	print 'known type', lknownType
	print 'each snp count', lsnpFreq
        print 'number of missing haplotype', nmiss

	#suppose the count of the new type is known as k
	for k in range(40):
		#for each possible number of count of the new type, solve integer linear programming
	        model = LpProblem('minor allele minimization problem', LpMinimize)
	
	        knownType_var = LpVariable.dicts('the index of each known haplotype', lknownType, lowBound = 1, upBound = 1000, cat = 'Integer')
	        newTypeSnp_var = LpVariable.dicts('each allele position of the new haplotype', range(nsnp),  cat = 'Binary')

                model += lpSum(newTypeSnp_var), 'Z'

		for j in range(nsnp):
			model += lpSum(knownType_var[stype] for stype in lknownType if stype[j]==lalt[j]) + k*newTypeSnp_var[j] == lsnpFreq[j]

		model += lpSum(knownType_var[i] for i in lknownType) + k == ntotal

		bbreak = False
		if model.solve() == 1:
			bbreak = True
			#print 'k', k
			#print model.solve()
			#print LpStatus
			for i in range(ntype):
				if i not in lknownType:
					continue
				if knownType_var[i].value() == 0.0:
					bbreak = False
			cnt1 = 0
			for j in range(nsnp):
				#print 'new type snp pos', j, newTypeSnp_var[j].value()
				if newTypeSnp_var[j].value() == None:
					bbreak = False
					continue
				elif newTypeSnp_var[j].value() >0 and newTypeSnp_var[j].value()<1:
					bbreak = False
					continue
				elif newTypeSnp_var[j].value()==1:
					cnt1 += 1
				else:	
					pass
			if cnt1 == 0:
				bbreak = False
		if bbreak:
			return k, knownType_var, newTypeSnp_var
	

	return 0,{},{}


def findNewType_2(lknownType, lsnpFreq, lalt, dtype, nmiss1, nmiss2):

        ltypes = dtype.keys()
        ntype = len(ltypes)
        print 'number of haplotypes', ntype+2
        for i in range(ntype):
                print 'real number of type', i, int(round(dtype.values()[i]*ntotal, 0))
        nsnp = len(lsnpFreq)
        print 'number of snps in the block, length of vector b', nsnp

        print 'known type', lknownType
        print 'each snp count', lsnpFreq


        #suppose the count of the new type is known as k
        for k1 in range(20):
                for k2 in range(20):
                #for each possible number of count of the new type, solve integer linear programming
                        model = LpProblem('minor allele minimization problem', LpMinimize)
        
                        knownType_var = LpVariable.dicts('the index of each known haplotype', lknownType, lowBound = 1, upBound = 1000, cat = 'Integer')
                        newTypeSnp_var_1 = LpVariable.dicts('each allele position of the new haplotype 1', range(nsnp),  cat = 'Binary')
                        newTypeSnp_var_2 = LpVariable.dicts('each allele position of the new haplotype 2', range(nsnp),  cat = 'Binary')
       

                        model += lpSum(newTypeSnp_var_1) + lpSum(newTypeSnp_var_2), 'Z'

                        for j in range(nsnp):
                                model += lpSum(knownType_var[stype] for stype in lknownType if stype[j]==lalt[j]) + k1*newTypeSnp_var_1[j] + k2*newTypeSnp_var_2[j] == lsnpFreq[j]

                        model += lpSum(knownType_var[i] for i in lknownType) + k1 + k2 == ntotal
        
                        bbreak = False
                        if model.solve() == 1:
                                bbreak = True
                                for i in range(ntype):
                                        if i not in lknownType:
                                                continue
                                        if knownType_var[i].value() == 0.0:
                                                bbreak = False
                                cnt1 = 0
                                cnt2 = 0
                                for j in range(nsnp):
                                        if newTypeSnp_var_1[j].value() == None or newTypeSnp_var_2[j].value() == None:
                                                bbreak = False
                                                continue
                                        elif newTypeSnp_var_1[j].value() >0 and newTypeSnp_var_2[j].value() >0 and newTypeSnp_var_1[j].value()<1 and newTypeSnp_var_2[j].value()<1:
                                                bbreak = False
                                                continue
                                        elif newTypeSnp_var_1[j].value()==1:
                                                cnt1 += 1
                                        elif newTypeSnp_var_2[j].value()==1:
                                                cnt2 += 1
                                        else:
                                                pass
                                if cnt1 == 0 or cnt2 == 0:
                                        bbreak = False
                        if bbreak:
                                return k1, k2, knownType_var, newTypeSnp_var_1, newTypeSnp_var_2
        
                        #print 'solution status', LpStatus[model.status]
                        #print 'model', model


        return 0,0,{},{},{}






if __name__ == '__main__':


	pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/sampling/newHaplotype/hapblock_leastFrequent_2Type.pk', 'rb')
	Lblock = pk.load(pkl_file)
	pkl_file.close()


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

	print 'ntotal', ntotal


	print 'total number of candidate blocks:', len(Lblock)

	cnt = 0
	lsolve = []
        for lblock in Lblock:
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
		lsnpFreq = [0.0 for i in range(len(lmarker))]
		lalt = []
		lref = []
                nsnp = len(lmarker)
		for i in range(len(lmarker)):
			idx = lmarker[i]
			idx = idx - 1
			sid = lrsidx[idx]
			lid = ['10'] + sid.split()[::-1]
			skey = ' '.join(lid)
			#print 'marker key', skey
			lskey.append(skey) #store the key to get the corresponding freq for each snp
			lallele = dsnpAllele[skey]
			#print 'lallele', lallele
			salt = dcode[lallele[-1]]
			lalt.append(salt)
			sref = dcode[lallele[0]]
			lref.append(sref)
			#lsnpFreq.append(dsnpFreq[skey]) #store the minor allele freq for each snp
	
			#sotre the minor allele freq for each snp for iterating each haplotype
			for stype in dtype.keys():
				#print 'sref, salt', sref, salt
				#print 'stype', stype[i]	
				if stype[i] == salt:	
					lsnpFreq[i] += dtype[stype] 				

		print 'dtype', dtype

		min_freq = min(dtype.values())
		print 'min freq', min_freq
		tidx_missing = dtype.values().index(min_freq)
                ltype_missing = []
                stype_missing = dtype.keys()[tidx_missing]
		ltype_missing.append(stype_missing)
                del dtype[stype_missing]
                min_freq_2 = min(dtype.values())
                print '2nd min freq',min_freq_2
                tidx_missing = dtype.values().index(min_freq_2)
                stype_missing_2 = dtype.keys()[tidx_missing]
                ltype_missing.append(stype_missing_2)
                del dtype[stype_missing_2]


		lsnpFreq = [int(round(xxx*ntotal,0)) for xxx in lsnpFreq] 

                k, knownType_var, newTypeSnp_var = findNewType(dtype.keys(), lsnpFreq, lalt, dtype, (int(round(min_freq*ntotal,0))+int(round(min_freq_2*ntotal,0)))) 

                #print 'k,knownTye_var, newTypeSnp_var', k, knownType_var, newTypeSnp_var


                if k!= 0:
                        print 'NEW TYPE FOUND!!!!!'
                        print '# of new type', k
                    
                        print 'false detection: new type error'
                        continue
        

                elif k == 0:
                        k1 = 0
                        k2 = 0
                        k1,k2, knownType_var, newTypeSnp_var_1, newTypeSnp_var_2  = findNewType_2(dtype.keys(), lsnpFreq, lalt, dtype, int(round(min_freq*ntotal,0)), int(round(min_freq_2*ntotal,0)))


                        #print 'k1,k2, knownType_var, newTypeSnp_var', k1,k2, knownType_var, newTypeSnp_var_1, newTypeSnp_var_2
		if k1 != 0 and k2 !=0:
			bt = True
			print 'NEW TYPE FOUND!!!!!'
			print '# of new type1', k1
                        print '# of new type2', k2
			snewType_1 = ''
                        for i in newTypeSnp_var_1:
                                print 'snp', i, 'value', newTypeSnp_var_1[i].value()
				if newTypeSnp_var_1[i].value() == 0.0:
					snewType_1 = snewType_1 + lref[i]
				elif newTypeSnp_var_1[i].value() == 1.0:
					snewType_1 = snewType_1 + lalt[i]
				else:
					print 'ERROR!!'
                        snewType_2 = ''
                        for i in newTypeSnp_var_2:
                                print 'snp', i, 'value', newTypeSnp_var_2[i].value()
                                if newTypeSnp_var_2[i].value() == 0.0:
                                        snewType_2 = snewType_2 + lref[i]
                                elif newTypeSnp_var_2[i].value() == 1.0:
                                        snewType_2 = snewType_2 + lalt[i]
                                else:
                                        print 'ERROR!!'
			print 'stype missing', stype_missing, stype_missing_2
                        print 'snewtype', snewType_1, snewType_2

                        output = open('solutionOuts_2type/block'+str(cnt)+'.pk', 'wb')
                        pk.dump([lblock, knownType_var, newTypeSnp_var_1, newTypeSnp_var_2, snewType_1, snewType_2, k1, k2], output)
                        output.close()


                        #print set([snewType_1.strip(), snewType_2.strip()]), set([stype_missing.strip(), stype_missing_2.strip()])
			if set([snewType_1.strip(), snewType_2.strip()]) != set([stype_missing.strip(), stype_missing_2.strip()]):
				bt = False
				print 'false detection: new type error'
                                print 'false detection: type count error'
                                continue
			for i in knownType_var:
				calp = knownType_var[i].value()
				realp = int(round(dtype[str(i)]*ntotal, 0))

				print 'knownType', i ,'count', knownType_var[i].value()
                		print 'real number of type', i, int(round(dtype[str(i)]*ntotal, 0))
				if abs(calp - realp) >5:
					bt = False
					print 'false detection: type count error'
                                        
                        if bt:
                            if abs(k1+k2-int(round(min_freq*ntotal,0))-int(round(min_freq_2*ntotal,0))) >10:
                                        bt = False
                                        print 'false detection: type count error'
                                        continue

    

			if bt:
				print 'success!!'
				#output = open('solutionOuts_2type_revise/block'+str(cnt)+'.pk', 'wb')
				#pk.dump([lblock, knownType_var, newTypeSnp_var_1, newTypeSnp_var_2], output)
				#output.close()
			
	
		#sys.exit()

	

