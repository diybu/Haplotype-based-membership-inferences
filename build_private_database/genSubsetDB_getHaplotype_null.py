
import sys
import os
import argparse
import random
import time
import pickle as pk
import copy
import numpy as np
import scipy.stats as st#import chisquare
import matplotlib.pyplot as plt
#import cPickle as pickle


def readfile(sfile):

	#lsubBlock = [[0 for j in range(len(Lskey))] for i in range(1000)]
	drsids = {}
        #print 'lmark', lmark
        with open(sfile, 'r') as f:
                for l in f:
                        l = l.strip()
                        if l.startswith('##'):
                                continue

                        if l.startswith('#CHROM'):
				r = l.split()
				N = len(r)
                                continue


                        r = l.split()
			#srsid = ' '.join(r[:3])
			srsid = r[2]			

			dtransSnp  = {'0': dcode[r[3]], '1': dcode[r[4]]}                               
			#ltypes = [] 
			stypes = ''
			for idx in range(9, len(r)):
				ldeploid = r[idx].strip().split('|')
				#ltypes.append(dtransSnp[ldeploid[0]])
				#ltypes.append(dtransSnp[ldeploid[1]])
				stypes = stypes + dtransSnp[ldeploid[0]]
				stypes = stypes + dtransSnp[ldeploid[1]]
			drsids[srsid] = stypes



	return drsids	

			



if __name__ == '__main__':

	dcode = {'A':'1', 'C':'2', 'G':'3', 'T':'4'}
	path = '/data2/diybu/beaconHaplopAttack/haploblock/'



        pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/samplePrivateDB/dsnpAllele_null.pk', 'rb')
        dsnpAllele = pk.load(pkl_file)
        pkl_file.close()






	dsnpTypes = readfile('chr3Subset.vcf')

        cntType = 0
        cntAll = 0
        cntNew = 0
        
	LblockSnp = []
	for filename in os.listdir('/data2/diybu/beaconHaplopAttack/haploblock/hapblockTarget_multiBlock_null'):

                pkl_file = open('/data2/diybu/beaconHaplopAttack/haploblock/hapblockTarget_multiBlock_null/'+filename, 'rb')
                lblock = pk.load(pkl_file)
                pkl_file.close()

                lmarker = lblock[0] #idx of marker, start from 1
                dtype = lblock[1] #type & freq A:1, C:2, G:3, T:4
                fname = lblock[-1].split('.')[0]

                #print dtype

                cntAll += len([xxx for xxx in dtype.values() if xxx < 0.1])
                
                #get rsid
                lrsidx = []
                with open('/data2/diybu/beaconHaplopAttack/haploblock/genPed_null/'+fname+'.info', 'r') as f:
                        for l in f:
                                sid = l.strip()
                                lrsidx.append(sid)

		
		#lskey = []
                lalt = []
		ltypes = ['' for xx in range(200)]
                for i in range(len(lmarker)):
			
                        idx = lmarker[i]
                        idx = idx - 1
                        sid = lrsidx[idx]
                        lid = ['3'] + sid.split()[::-1]
                        skey = ' '.join(lid)
			srsid = lid[-1]
                        #lskey.append(skey) #store the key to get the corresponding freq for each snp
			print 'srsid, sref, salt', srsid, dcode[dsnpAllele[skey][0]], dcode[dsnpAllele[skey][1]]
			#print 'dsnpTypes', dsnpTypes

			for j in range(200):
				ltypes[j] = ltypes[j] + dsnpTypes[srsid][j]


                ltypes = list(set(dtype.keys())-set(ltypes))
                #cntType += len(ltypes)


                ltarget = []
                for stype in ltypes:
                        if stype not in dtype or dtype[stype] >= 0.1:
                                continue
                        else:
                                print 'dtype[stype]',dtype[stype]
                                #if stype not in ltarget and dtype[stype]<0.1:
                                ltarget.append(stype)
                                cntType += 1
                
                if ltarget:
                        sout = filename.split('.')[0]+'_nullTypes.pk'
        	        output = open('nullType_for_nullDist/'+sout, 'wb')
                        pk.dump([lblock, ltarget], output)
                        output.close()

        print '# of target types w/ prior prob < 0.1 in the DB:', cntType
        print '# of target types w/ prior prob < 0.1 in public resource', cntAll
