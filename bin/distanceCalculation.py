import os, sys
import glob

### Supporting script for maximum likelihood calculations ###

def calculateProteinDistances(species1,species2,config):
    root_dir = config.path_distance_work_dir
    omaSeqs = config.path_oma_seqs
    omaPairs = config.path_oma_pairs
    concatAlignment = config.concat_alignments_scripts
    omaProteomesDir = ".."
    
    if not os.path.exists('log.txt'):
    	log = open('log.txt', 'w')
    else:
    	log = open('log.txt', 'a')
    
    workDir = species1+'_'+species2
    
    if not os.path.exists(workDir):
    	os.mkdir(workDir)
    
    # The directory which contains the fasta files
    faDir = os.path.abspath(workDir) + '/faDir'
    # The directory which contains the alignments
    alnDir = os.path.abspath(workDir) + '/alnDir'
    
    if not os.path.exists(faDir):
    	os.mkdir(faDir)
    
    if not os.path.exists(alnDir):
    	os.mkdir(alnDir)
    
    def preprocess():
    	print 'Preprocessing:\tParsing sequence pairs and aligning them..'
    	flag = True
    	c = 1
    	for line in open(omaPairs):
    		if species1 in line and species2 in line:
    			flag = False
    			print 'Sequence Nr.:', c
    			if species1 in line.split('\t')[0]:
    				species1Id = line.split('\t')[0]
    				species2Id = line.split('\n')[0].split('\t')[1]
    			else:
    				species2Id = line.split('\t')[0]
    				species1Id = line.split('\n')[0].split('\t')[1]
    
    
    			filename = faDir + '/seq_'+str(c)+'.fa'
    			fnew = open(filename, 'w')
    
    			if os.path.exists(omaProteomesDir + '/proteome_' + species1):
    				with open(omaProteomesDir + '/proteome_' + species1) as prot1:
    					for protLine in prot1:
    						if '>' in protLine and species1Id in protLine:
    							fnew.write(protLine[:6] + '\n' + prot1.next().replace('*', ''))
    							break
    
    			else:
    				print 'Searching in file: ', omaSeqs, species1
    				with open(omaSeqs) as omaseq:
    					for line2 in omaseq:
    						if '>' in line2 and species1Id in line2:
    							fnew.write(line2[:6] + '\n' + omaseq.next().replace('*', ''))
    							break
    
    			if os.path.exists(omaProteomesDir + '/proteome_' + species2):
    				with open(omaProteomesDir + '/proteome_' + species2) as prot2:
    					for protLine in prot2:
    						if '>' in protLine and species2Id in protLine:
    							fnew.write(protLine[:6] + '\n' + prot2.next().replace('*', ''))
    							break
    			else:
    				print 'Searching in file: ', omaSeqs, species2
    				with open(omaSeqs) as omaseq:
    					for line2 in omaseq:
    						if '>' in line2 and species2Id in line2:
    							fnew.write(line2[:6] + '\n' + omaseq.next().replace('*', ''))
    							break
    			fnew.close()
    			c += 1
    
    			os.system('muscle -quiet -in %s -out %s' %(filename, filename.replace('.fa', '.aln').replace(faDir, alnDir)))
    			#os.system('linsi %s > %s' %(filename, filename.replace('.fa', '.aln').replace(faDir, alnDir)))
    
    	if flag:
    		print 'WARNING! ERROR! Species missing!'
    		log.write('WARNING: Either %s or %s is missing in the OMA database' %(species1, species2))
    
    def postprocess():
    	print 'Preprocessing complete..\nConcating the alignments..'
    	concatFile = species1+'_'+species2+'.aln'
    	os.system('perl %s -in=%s -out=%s' %(concatAlignment, alnDir, concatFile))
    	print 'Postprocessing concatenated alignment file..'
    	postProcessConcatFile = concatFile.replace('.aln', '.dup.aln')
    	fnew = open(postProcessConcatFile, 'w')
    	with open(concatFile) as con:
    		for line in con:
    			if '>' in line:
    				seqId = line
    				sequence = con.next()
    				fnew.write(seqId + sequence + seqId.replace('>', '>dup_') + sequence)
    	fnew.close()
    	print 'Converting alignment file to phylip format..'
    	phyFile = postProcessConcatFile.replace('.aln', '.phy')
    	os.system('clustalw -convert -output=phylip -infile=%s -outfile=%s' %(postProcessConcatFile, phyFile))
    
    	print 'Performing likelihood mapping..'
    	puzzleParams = 'temp_puzzleParams.txt'
    	fnew = open(puzzleParams, 'w')
    	fnew.write(phyFile + '\nb\ne\nm\nm\nm\nm\nm\nm\ny\n')
    	fnew.close()
    	os.system('puzzle < %s' %puzzleParams)
    
    	print 'Parsing outfile and saving likelihood distance to'
    	result = species1+'_'+species2+'.lik'
    	fnew = open(result, 'w')
    	outdist = open(phyFile + '.dist').read().split('\n')
    	fnew.write(outdist[3].split()[1])
    	fnew.close()
    
    #	if os.path.exists(result) and not len(open(result).read().split('\n')) == 0:
    #		os.system('rm -rf %s' %faDir)
    #		os.system('rm -rf %s' %alnDir)
    #		os.system('rm %s*' %phyFile)
    #		os.system('rm -rf *.txt')
    #		os.system('rm -rf *.aln')
    
    os.chdir(workDir)
    preprocess()
    postprocess()
    os.chdir(rootDir)
