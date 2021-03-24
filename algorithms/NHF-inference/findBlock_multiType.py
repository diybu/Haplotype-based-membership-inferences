
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
import cPickle as pickle


def readfile(path):
        lfName = []
	global LhapblockTarget
	LhapblockTarget = []
	lpars = [] # store all length and freq for 15 haploblock candidates
	lpar = [1.0, 0] # store max length and freq for 15 haploblock candidates
        for root, dirs, files in os.walk(path):
                for fname in files:
                        pfname = path+fname
			#fname = '/data2/diybu/beaconHaplopAttack/haploblock/hapOut/blocks/10_9900000-10070000.GABRIELblocks'

			nstart = int(fname.split('-')[0].split('_')[1])
			nend = int(fname.split('-')[1].split('.')[0])
			#print nstart
			#continue

			with open(pfname, 'r') as f:
        			freq = 1.0
				Ltypes = []
				dtype = {}
				for l in f:
					if l.startswith('Multiallelic'):
						continue
					if l.startswith('BLOCK'):
						if dtype:
                                                        Ltypes.append([lmark, dtype, fname])
						r = l.strip().split('MARKERS:')[1].strip().split()
						lmark = map(int, r)
						#lmark = map(lambda x:x+nstart, lint)
						#if lmark[0] < nstart or lmark[-1] > nend:
							#print fname
							#print lmark
							#sys.exit()
						#lmark = r
						dtype = {}
						continue
					r = l.strip().split('|')[0].strip(')\t').split('(')
					k = r[0]
					v = float(r[1])
					dtype[k] = v
				Ltypes.append([lmark, dtype, fname])




###################### least frequent ###################################################################################
                        Lselect, freq = lowDomimant(Ltypes, freq)
			#print 'freq', freq
			if freq > 0.02:
				continue
                        if Lselect and len(LhapblockTarget) < 10000:
                                #print lpars
                                #print 'lpar', lpar
                                LhapblockTarget.append(Lselect)
                                lpars.append([freq, len(lpars)])
                                if lpar == [1.0, 0]:
                                        lpar = [freq, len(lpars)-1]
                                elif freq < lpar[0]:
                                        lpar = [freq, len(lpars)-1]
                        elif Lselect:
                                #print lpars
                                
                                if freq < lpar[0]:
					#print 'freq, lpar', freq, lpar
                                        index = lpar[1]
                                        LhapblockTarget[index] = Lselect
                                        lpar = [freq, index]
                                        lpars[index] = lpar


	print '# of select blocks', len(LhapblockTarget)
        fhapblock = open('hapblock_leastFrequent_3Type.pk', 'wb')
        pickle.dump(LhapblockTarget, fhapblock)
        fhapblock.close()

	return




	





def lowDomimant(Ltypes, freq):
#find the most robust haploblocks with rare haplotypes
        Lselect = []

        for lhapblock in Ltypes:
                #print lhapblock
                dhapblock = lhapblock[1]
                freqtmp = min(dhapblock.values())
                nmin = len([x for x in dhapblock.values() if x <=0.02])
                if freqtmp < freq and nmin >=3:
                        Lselect = [lhapblock[0], dhapblock, lhapblock[2]]
                        freq = freqtmp

        return Lselect, freq





if __name__ == '__main__':


	path = '/data2/diybu/beaconHaplopAttack/haploblock/hapOut/blocks/'
	dtype = readfile(path)
	#print dtype

        pass
