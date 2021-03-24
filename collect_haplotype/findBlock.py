
import sys
import os
import argparse
import random
import time
import pickle
import copy
import numpy as np
import scipy.stats as st#import chisquare
import math


def readfile(path):
        min_marker = float('inf')
        max_marker = -float('inf')
        cnt_target = 0
        test_num = 0 # record num of block candidates
        for root, dirs, files in os.walk(path):
                for fname in files:
			if test_num > 1000:
				break
                        pfname = path+fname

			nstart = int(fname.split('-')[0].split('_')[1])
			nend = int(fname.split('-')[1].split('.')[0])
			#print nstart
			#continue

			with open(pfname, 'r') as f:
				length = 0
        			freq = 0.0
				Ltypes = []
				dtype = {}
				for l in f:
					if l.startswith('Multiallelic'):
						continue
					if l.startswith('BLOCK'):
						if dtype:
                                                        Ltypes.append([lmark, ltypes, dtype, Llinkp, fname])
						r = l.strip().split('MARKERS:')[1].strip().split()
						lmark = map(int, r)
						#lmark = map(lambda x:x+nstart, lint)
						#if lmark[0] < nstart or lmark[-1] > nend:
							#print fname
							#print lmark
							#sys.exit()
						#lmark = r
						dtype = {}
						Llinkp = []
						ltypes = []
						continue

					r = l.strip().split('|')[0].strip(')\t').split('(')
					k = r[0]
					v = float(r[1])
					dtype[k] = v
					ltypes.append(k)
				
					if '|' not in l:
						continue
					llink = l.strip().split('|')[1].strip('|').split()
					Llinkp.append(llink)
					#print 'fname, Llinkp', fname, Llinkp
				Ltypes.append([lmark, ltypes, dtype, Llinkp, fname])


                        dtype = {}
                        Llinkp = []
                        lmark = []
                        ltypes = []
			print 'fname', fname
                        print '# of block in this file', len(Ltypes)
######extend haplotype until the prior prob smaller than threshold 0.1
			pprior = 1.0
			cnt = 0
                        cnt_begin = 0
			nmarker = 0		
			#print 'Ltypes', Ltypes
			
			dtype_prior = {}
			dtype_freq = {}
			dtype_idx = {}

                        while cnt < len(Ltypes)-1:
                                bnew = True
                                pprior = 1.0
                                nmarker = 0
                                dtype_prior = {}
                                dtype_freq = {}
                                dtype_idx = {}
                                print 'first while per file','cnt', cnt
        			while nmarker < 50 and pprior > 0.1 and cnt < len(Ltypes)-1:
        				print 'second while, per target block', 'cnt', cnt
        		
                                        
        				Lselect = Ltypes[cnt]
        				ltypes = Lselect[1]

                                        cnt += 1
        			
        				if len(ltypes) > 10:
        					break
        	
        				dtype = Lselect[2]
        				Llinkp = Lselect[3]
        				if Llinkp == []:
        					break
        				elif len(Llinkp[0]) > 10:
        					break
        
        				#print 'ltypes', ltypes
        				#print 'Llinkp', Llinkp
        				
        				if bnew:
                                                cnt_begin = cnt - 1
        					#print 'cnt', cnt
                                                #print 'dtype_prior', dtype_prior
        					for ii in range(len(ltypes)):		
        						stype = ltypes[ii]
        						freq = dtype[stype]
        						prior = (1-(1-freq)**1000)
                                                        #print 'prior', prior, freq
        						dtype_idx[str(ii)] = stype.strip()
        						dtype_prior[str(ii)] = math.log(prior)
        						dtype_freq[str(ii)] = math.log(freq)
        						for ij in range(len(Llinkp[ii])):			
        							plink = float(Llinkp[ii][ij])
        							if plink == 0.0:
        								continue
        							dtype_prior[str(ii)+str(ij)] = dtype_prior[str(ii)] + math.log(plink)
        							dtype_freq[str(ii)+str(ij)] = dtype_freq[str(ii)] + math.log(plink)
        							if dtype_prior[str(ii)+str(ij)] != 0.0 and pprior > dtype_prior[str(ii)+str(ij)]:
        								pprior = math.exp(dtype_prior[str(ii)+str(ij)])
                                                                        print 'pprior', pprior
        						del dtype_freq[str(ii)]
                                                	del dtype_prior[str(ii)]
                                                bnew = False
        					#print 'dtype_prior', dtype_prior
        				else:
        					#print 'cnt', cnt
        					#print 'dtype_prior', dtype_prior
        					for key in dtype_prior.keys():
        						prev = key[:-1]
        						ii = int(key[-1])
        						stype = ltypes[ii]
                                                        freq = dtype[stype]
                                                        prior = (1-(1-freq)**1000)
                                                        #print 'prior', prior, freq
        
        						dtype_idx[key] = dtype_idx[prev]+stype.strip()
                                                	dtype_prior[key] = dtype_prior[key]+math.log(prior)
                                                	dtype_freq[key] = dtype_freq[key]+math.log(freq)
                                                	for ij in range(len(Llinkp[ii])):
                                                	        plink = float(Llinkp[ii][ij])
        							#print 'plink', plink
        							if plink == 0.0:
        								continue
                                                	        dtype_prior[key+str(ij)] = dtype_prior[key] + math.log(plink)
                                                	        dtype_freq[key+str(ij)] = dtype_freq[key] + math.log(plink)
                                                	        if dtype_prior[key+str(ij)] != 0.0 and pprior > dtype_prior[key+str(ij)]:
                                                	                pprior = math.exp(dtype_prior[key+str(ij)])
                                                                        print 'pprior', pprior
        						del dtype_prior[key]
                                                	del dtype_freq[key]
        					
        				nmarker += len(stype.strip())
        				#cnt += 1	
        
                                        print 'nmarker', nmarker
        				print 'pprior', pprior
        				print 'cnt', cnt
        		
        			if pprior > 0.1:
        				continue
        
        			if pprior == 0.0:
        				continue
        
        			if nmarker < 50:
        				continue 
        
        			
        			if not bnew: #add linked types prob to freq & pprior
        				Lselect = Ltypes[cnt]
                                        ltypes = Lselect[1]                                
                                        dtype = Lselect[2]
        				#print 'final cnt', cnt
        				#print 'ltypes', ltypes
        			
        				#print 'dtype_prior', dtype_prior
        	
        				for key in dtype_prior.keys():
        					ii = int(key[-1])
        					stype = ltypes[ii]
                                                freq = dtype[stype]
                                                prior = (1-(1-freq)**1000)
        
        					prev = key[:-1]
        					dtype_idx[key] = dtype_idx[prev]+stype.strip()
                                        	dtype_prior[key] = dtype_prior[key]+math.log(prior)
                                                dtype_freq[key] = dtype_freq[key]+math.log(freq)
        				Lmarker = []
                                        for xx in range(cnt_begin, cnt+1):
                                                Lselect = Ltypes[xx]
                                                Lmarker += Lselect[0]
                                else:
        			        Lmarker = []
        			        for xx in range(cnt_begin, cnt):
        				        Lselect = Ltypes[xx]
        				        Lmarker += Lselect[0]
        
        		        Lselect = []	
                                Ltypes = []
        
        			dtype_prior_clean = {}
        			dtype_freq_clean = {}
        			#print 'dtype_prior', dtype_prior
        			for key in dtype_idx.keys():
        				if key not in dtype_prior:
        					continue
        				print 'key', key, len(dtype_idx[key]), len(Lmarker)
        				skey = dtype_idx[key]
        				if dtype_prior[key] >= 0.0 or dtype_freq[key] >= 0.0:
        					continue
        				dtype_prior_clean[skey] = math.exp(dtype_prior[key])
                                        if dtype_prior_clean[skey] < 0.1:
                                                cnt_target += 1
        				dtype_freq_clean[skey] = math.exp(dtype_freq[key])
       

                                if not (dtype_prior) or (not dtype_freq):
                                        continue
                                print 'dtype_prior', dtype_prior_clean, dtype_prior
                                print 'dtype_freq', dtype_freq_clean, dtype_freq

                                dtype_prior = {}
                                dtype_freq = {}
                                dtype_idx = {}
        
                                nmarker = len(Lmarker)
                                print 'nmarker', nmarker
                                print 'pprior', pprior
                                if pprior < 0:
                                        sys.exit()
                                if nmarker < min_marker:
                                        min_marker = nmarker
                                if nmarker > max_marker:
                                        max_marker = nmarker                                
        
        			Lcandidate = [Lmarker, dtype_freq_clean, dtype_prior_clean, fname]			
        
                                Lmarker = []
                                dtype_prior_clean = {}
                                dtype_freq_clean = {}
        
        			test_num += 1
        		
        
        			#LhapblockTarget.append(Lcandidate)
        
                                fhapblock = open('hapblockTarget_multiBlock/block'+str(test_num)+'.pk', 'wb')
                                pickle.dump(Lcandidate, fhapblock)
                                fhapblock.close()
        	
                                Lcandidate = []
        
        print '# of target rare types', cnt_target

        print 'min marker', min_marker
        print 'max marker', max_marker
	return












if __name__ == '__main__':


	path = '/data2/diybu/beaconHaplopAttack/haploblock/hapOut/blocks/'
	dtype = readfile(path)
	#print dtype

        pass
