#
# Module to perform preprocessing steps in protTrace workflow
#
import os, sys
import glob
import configure
import hamstr_search
import treeReconstruction
import transformAlignment
import dendropy
import subprocess
import random
import time
from operator import itemgetter
import xml.etree.ElementTree as ET
from optparse import OptionParser
from multiprocessing import Pool

def Preprocessing(prot_id, querySeq, config_file):
	# Store the current working directory into a variable
	rootDir = os.getcwd()

	print 'Prot_id: ', prot_id

	# Creating instance of class configure
	# Saves the information provided by the program configuration file
	prot_config = configure.setParams(config_file)

	# Declares global variables which will be used by all methods in the module
	global include_paralogs, hamstr_env, cache, cache_dir, work_dir, omaIdFile, orth_file, fas_file, aln_file, phy_file, id_file, proteome_file, tree_file, trans_file, hmm_file, xml_file, REvolver_output_dir, species_id, indel_file, scale_file, ortholog_tree_reconstruction, nr_processors, hamstr_oma_map_file, delTemp, blastp

	# Getting information from the configuration file
	# Setting the names for the protTrace temporary and output files	
	species_id = prot_config.species
	proteome_file = 'proteome_' + prot_id
	id_file = 'ogIds_' + prot_id + '.txt'
	work_dir = prot_config.path_work_dir + '/' + prot_id
	cache_dir = prot_config.path_cache
	omaIdFile = work_dir + '/omaId.txt'
	orth_file = 'ogSeqs_' + prot_id + '.fa'
	fas_file = 'ogSeqs_' + prot_id + '.fasScore'
	domain_archi_file = 'ogSeqs_' + prot_id + '.domains'
	aln_file = 'ogSeqs_' + prot_id + '.aln'
	phy_file = 'ogSeqs_' + prot_id + '.phy'
	tree_file = 'RAxML_bestTree.' + prot_id 	###	CHANGE HERE IF RAxML OUTPUT NAME CHANGES 	###
	trans_file = 'ogSeqs_' + prot_id + '.trans'
	hmm_file = prot_id + '.hmm'
	xml_file = 'revolver_config_' + prot_id + '.xml'
	REvolver_output_dir = work_dir + '/REvolver_output/'
	indel_file = 'indel_' + prot_id
	scale_file = 'scale_' + prot_id
	delTemp = prot_config.delete_temp
	cache = prot_config.reuse_cache
	hamstr_env = prot_config.hamstr_environment
	include_paralogs = prot_config.includeParalogs
	ortholog_tree_reconstruction = prot_config.phylogeneticTreeReconstruction
	nr_processors = prot_config.nr_processors
	hamstr_oma_map_file = prot_config.hamstr_oma_tree_map
	blastp = prot_config.blastp

	# Creates a working directory where temporary and output files will be stored
	if not os.path.exists(work_dir):
		print '##### Creating working directory:\n', work_dir
		try:
			os.mkdir(work_dir)
		except:
			sys.exit('ERROR: Working directory cannot be created!')

	# Change current working directory
	os.chdir(work_dir)

	# Parse proteome of the input species given in program configuration file
	# The proteome is extracted from the OMA database sequences file
	startProcessTime = time.time()
	parseProteome(species_id, prot_config.path_oma_seqs, prot_config.makeblastdb, proteome_file, prot_config.hamstr_oma_tree_map, prot_config.hamstr, hamstr_env, prot_config.search_oma_database)
	proteome_file = os.path.abspath(proteome_file)
	print '#####\tTIME TAKEN: %s mins\tSpecies %s gene set preparation#####' %((time.time() - startProcessTime) / 60, species_id) 	

	# Search ortholog groups for the input OMA id
	# In case of fasta sequences, first the OMA id is parsed from the OMA database
	# followed with the extraction of the OMA group
	if prot_config.orthologs_prediction:

		if cache and os.path.exists(orth_file):
			print 'Orthologs file exist. Reusing it.'
		elif cache and os.path.exists(cache_dir + '/' + orth_file):
			print 'Orthologs file found in cache_dir. Reusing it.'
			os.system('cp %s %s' %(cache_dir + '/' + orth_file, orth_file))
		else:
			# Check if OMA database has to be accessed for an ortholog search.
			if prot_config.search_oma_database:
				startProcessTime = time.time()
				# Find OMA orthologs groups if any	
				run = findOmaOrthologs(prot_id, querySeq, prot_config.path_oma_group, prot_config.path_oma_seqs, proteome_file, prot_config.formatdb, prot_config.blastp, delTemp, species_id)		

				# Search for the ortholog sequences for the respective OMA orthologs group
				# For all the OMA ids in the OMA group, extract sequences from OMA database sequences file
				if run == 2:
					findOmaSequences(prot_id, prot_config.path_oma_seqs, species_id, prot_config.hamstr_oma_tree_map)
					print '#####\tTIME TAKEN: %s mins\tOrthologs search in OMA database.\t#####' %((time.time() - startProcessTime) / 60)
				else:
					print '##### Preparing ortholog file #####'
					fOrth = open(orth_file, 'w')
					fOrth.write('>' + species_id + '\n' + querySeq)
					fOrth.close()

			# HaMStR / HaMStR-OneSeq search
			# The orthologs sequences by OMA is used as core-ortholog set
			try:
				f = open(orth_file).read().split('\n')
				# Run HaMStR search if 2 or more sequences are present. Otherwise, run HaMStROneSeq search if only 1 sequence is present
				if len(f) > 3:
					if prot_config.run_hamstr:
						startProcessTime = time.time()
						success = hamstr_search.main(prot_config.hamstr, orth_file, prot_id, prot_config.hamstr_oma_tree_map, prot_config.formatdb, prot_config.blastp, delTemp, hamstr_env, include_paralogs, nr_processors)
						if success:
							print '#####\tTIME TAKEN: %s mins\tHaMStR search in local genome directory.\t#####' %((time.time() - startProcessTime) / 60)
							os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))
						else:
							if prot_config.run_hamstrOneSeq:
								print '##### HaMStROneSeq search for orthologs #####'
								startProcessTime = time.time()

								# Read the orthologs file and limit it to just the query species id and sequence
								ortholog_temp = open(orth_file).read().split('\n')
								rewrite_orth_file = open(orth_file, 'w')
								for orthLines in range(len(ortholog_temp) - 1):
									if '>' in ortholog_temp[orthLines] and species_id in ortholog_temp[orthLines]:
										rewrite_orth_file.write(ortholog_temp[orthLines] + '\n' + ortholog_temp[orthLines + 1])
										break
								rewrite_orth_file.close()
								run_hamstrOneSeq(prot_config.hamstr, os.path.abspath(orth_file), prot_config.hamstr_oma_tree_map, prot_id, prot_config.formatdb, prot_config.blastp, proteome_file, delTemp, hamstr_env, include_paralogs)
								print '#####\tTIME TAKEN: %s mins\tHaMStR-OneSeq#####' %((time.time() - startProcessTime) / 60)
								os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))

					elif prot_config.run_hamstrOneSeq:
						startProcessTime = time.time()
						print '##### HaMStROneSeq search for orthologs #####'
						
						# Read the orthologs file and limit it to just the query species id and sequence
						ortholog_temp = open(orth_file).read().split('\n')
						rewrite_orth_file = open(orth_file, 'w')
						inputTaxaSet = open('inputTaxaSet_oneSeq.txt', 'w')
						for orthLines in range(len(ortholog_temp) - 1):
							if '>' in ortholog_temp[orthLines] and species_id in ortholog_temp[orthLines]:
								rewrite_orth_file.write(ortholog_temp[orthLines] + '\n' + ortholog_temp[orthLines + 1])
							elif '>' in ortholog_temp[orthLines] and not species_id in ortholog_temp[orthLines]:
								inOmaId = ortholog_temp[orthLines].split()[0][1:]
								for mapLine in open(prot_config.hamstr_oma_tree_map):
									if inOmaId in mapLine:
										inputTaxaSet.write(mapLine.split()[0] + '\n')
										break

						inputTaxaSet.close()
						rewrite_orth_file.close()
						run_hamstrOneSeq(prot_config.hamstr, os.path.abspath(orth_file), prot_config.hamstr_oma_tree_map, prot_id, prot_config.formatdb, prot_config.blastp, proteome_file, delTemp, hamstr_env, include_paralogs)
						print '#####\tTIME TAKEN: %s mins\tHaMStR-OneSeq#####' %((time.time() - startProcessTime) / 60)
						os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))
				

				elif len(f) > 0 and len(f) < 4:
					if prot_config.run_hamstrOneSeq:
						print '##### HaMStROneSeq search for orthologs #####'
						startProcessTime = time.time()
						if cache and os.path.exists(cache_dir + '/' + orth_file):
							print 'Pre-HaMStR computed orthologs file found in Cache for re-use!'
							os.system('cp %s %s' %(cache_dir + '/' + orth_file, orth_file))
						else:
							run_hamstrOneSeq(prot_config.hamstr, os.path.abspath(orth_file), prot_config.hamstr_oma_tree_map, prot_id, prot_config.formatdb, prot_config.blastp, proteome_file, delTemp, hamstr_env, include_paralogs)
							print '#####\tTIME TAKEN: %s mins\tHaMStR-OneSeq#####' %((time.time() - startProcessTime) / 60)
							os.system('cp %s %s' %(orth_file, cache_dir + '/' + orth_file))		

				else:
					sys.exit('ERROR: No sequence found in OMA sequences! The ortholog sequences file is empty!')
			except IOError:
				sys.exit('ERROR: Orthologs sequences file is invalid!')

			except KeyboardInterrupt:
				sys.exit('Keyboard interruption by user!!!')

	# Calculate FAS scores for identified orthologs
	if prot_config.fas_score:
		if cache and os.path.exists(cache_dir + '/' + fas_file):
			print 'Pre-computed FAS score file found in cache. Skipping FAS calculations.'
			os.system('cp %s %s' %(cache_dir + '/' + fas_file, fas_file))
		else:
			print '##### Calculating FAS scores for identified orthologs #####'
			startProcessTime = time.time()
			calculateFAS(work_dir, prot_config.hamstr, prot_config.fas_annotations, os.path.abspath(orth_file), os.path.abspath(fas_file), os.path.abspath(domain_archi_file), prot_config.hamstr_oma_tree_map, prot_config.blastp, prot_id, species_id, delTemp, hamstr_env)
			print '#####\tTIME TAKEN: %s mins\tFAS#####' %((time.time() - startProcessTime) / 60)
			os.system('cp %s %s' %(fas_file, cache_dir))

	# Performs MSA on the orthologs sequences
	if prot_config.perform_msa:
		print '##### Performing MSA of the orthologs sequences #####'
		startProcessTime = time.time()
		performMSA(prot_config.msa, prot_config.clustalw)
		print '#####\tTIME TAKEN: %s mins\tMAFFT#####' %((time.time() - startProcessTime) / 60)

	# Calls tree reconstruction module which generates tree using degapped alignment
	# and also calculates the scaling factor based on maximum likelihood distance between species
	if prot_config.calculate_scaling_factor:
		if cache and os.path.exists(scale_file):
			print 'Pre-computed scaling factor found for re-use!'
		else:
			print '##### Tree reconstruction and scaling factor calculation #####'
			startProcessTime = time.time()
			treeReconstruction.main(prot_config.tree_reconstruction, prot_config.msa, prot_config.clustalw, prot_config.degapping, orth_file, prot_config.aa_substitution_matrix, prot_id, prot_config.treePuzzle, prot_config.parameters_treePuzzle, prot_config.hamstr_oma_tree_map, prot_config.species_MaxLikMatrix, scale_file, tree_file, delTemp, prot_config.default_scaling_factor, cache_dir, ortholog_tree_reconstruction, nr_processors)
			print '#####\tTIME TAKEN: %s mins\tRAxML#####' %((time.time() - startProcessTime) / 60)

	# Calculate indels
	if prot_config.calculate_indel:
		if cache and os.path.exists(indel_file):
			print 'Pre-computed indel found for re-use!'
		else:

			# Transform alignment
			print '##### Transforming MSA based on indel blocks #####'
			alignmentLength = 0
			try:
				alignmentLength = transformAlignment.main(phy_file, trans_file)
			except:
				pass
			calculateIndels(tree_file, trans_file, alignmentLength, prot_config.iqtree24, prot_config.default_indel, prot_config.default_indel_distribution)

	# Domain constraint file for REvolver
	if prot_config.traceability_calculation:
		# Creates a output directory for REvolver
		if not os.path.exists(REvolver_output_dir):
			print '##### Creating REvolver output directory:\n', REvolver_output_dir
			try:
				os.mkdir(REvolver_output_dir)
			except:
				sys.exit('ERROR: REvolver output directory cannot be created!')

		print '##### Generating domain constraints for REvolver #####'
		hmmscan(prot_config.hmmscan, orth_file, prot_config.pfam_database, hmm_file, prot_id, species_id)

		# Prepare XML config file to be used as an input for REvolver
		print '##### Preparing XML configuration file for REvolver #####'
		if os.path.exists(scale_file):
			f = open(scale_file).read().split('\n')
			scaling_factor = f[0]
		else:
			print 'WARNING: Scaling factor file not found. Using default value:', prot_config.default_scaling_factor
			scaling_factor = prot_config.default_scaling_factor
			# Writing randomly generated indels into the file (Randome pick between 0.50 to 2.00)
			writeScale = open(scale_file, 'w')
			writeScale.write(scaling_factor)
			writeScale.close()

		if os.path.exists(indel_file):
			f = open(indel_file).read().split('\n')
			indel = f[0]
			p = f[1]
		else:
			print 'WARNING: Indel file not found. Using default value:', prot_config.default_indel
			indel = prot_config.default_indel
			p = prot_config.default_indel_distribution
			# Writing randomly generated indels into the file (Random pick for indel between 0.05 to 1.9; for distribution between 0.0 to 0.9)
			# The random value limit has been set as observed from the large scale indels distribution
			writeIndel = open(indel_file, 'w')
			writeIndel.write(indel + '\n' + p)
			writeIndel.close()

		prepareXML(xml_file, prot_config.pfam_database, prot_config.hmmfetch, prot_config.aa_substitution_matrix, indel, p, scaling_factor, prot_config.simulation_tree, prot_id, hmm_file, REvolver_output_dir)
	
	os.chdir(rootDir)

# Feature Architecture parser
def parseFAS(fasOutputFile,prot_id, sp1, sp2):
	domainArchiData = []
	name = prot_id
	#if not os.path.exists(domainArchiFile):
	#	out = open(domainArchiFile, "w")
	#else:
	#	out = open(domainArchiFile, "a")
	xmltree = ET.parse(fasOutputFile)
	root = xmltree.getroot()
	set_dic = {}
	if root.attrib["direction"] == "single-->set":
		
		for protein in root:
			if protein.tag == "single_protein":
				infodump = (sp1 + ':' + protein.attrib["id"], 0, [])
		for protein in root:
			if protein.tag == "set_protein":
				pid = sp2 + ':' + protein.attrib["id"]
				p_features = []
				path = []
				for arc in protein:
					if arc.tag == "architecture":
						for feature in arc:
							ftype = feature.attrib["type"]
							weight = feature.attrib["weight"]
							for instance in feature:
								#out.write(name + "#" + pid + "#" + infodump[0] + "\t" + pid + "\t" + str(ftype) + "\t" + str(instance.attrib["start"]) + "\t" + str(instance.attrib["end"]) + "\t" + str(weight) + "\n")
								domainArchiData.append(name + "#" + pid + "#" + infodump[0] + "\t" + pid + "\t" + str(ftype) + "\t" + str(instance.attrib["start"]) + "\t" + str(instance.attrib["end"]) + "\t" + str(weight) + "\n")
					elif arc.tag == "path":
						for feature in arc:
							ftype = feature.attrib["type"]
							weight = feature.attrib["corrected_weight"]
							for instance in feature:
								#out.write(name + "#" + pid + "#" + infodump[0] + "\t" + infodump[0] + "\t" + str(ftype) + "\t" + str(instance.attrib["start"]) + "\t" + str(instance.attrib["end"]) + "\t" + str(weight) + "\n")
								domainArchiData.append(name + "#" + pid + "#" + infodump[0] + "\t" + infodump[0] + "\t" + str(ftype) + "\t" + str(instance.attrib["start"]) + "\t" + str(instance.attrib["end"]) + "\t" + str(weight) + "\n")
	#out.close()
	return domainArchiData

### FAS Annotations computation - annotation.pl script works here ###
def retrieve_FAS_annotations(fasta):
	# Read the orthologs file and create the FAS annotations for individual protein sequences	
	oma_ids = []
	mapDict = {}
	fasCacheData = []
	logDict = {}

	omaId = fasta.split('\n')[0][1:]

	if omaId.split('_')[1] == species_id:
		oma_ids.append(omaId.split('_')[1])
	else:
		oma_ids.append(omaId)
	sequence = fasta.split('\n')[1]

	for mapLine in open(hamstr_oma_map_file):
		if omaId.split('_')[1] == mapLine.split()[-1]:
			hamstrId = mapLine.split()[0]
			if not omaId in mapDict.keys():
				if omaId.split('_')[1] == species_id:
					mapDict[omaId.split('_')[1]] = hamstrId
				else:
					mapDict[omaId] = hamstrId
			break

	species_anno_dir = speciesAnnoDir + '/' + hamstrId
	if not os.path.exists(species_anno_dir):
		print 'ERROR: Species %s annotation files not found. Please check proper path is provided in configuration file.' %hamstrId
		sys.exit()

	if omaId.split('_')[1] == species_id:
		protein_temp_anno_dir = temp_fasAnno + '/' + omaId.split('_')[1]
	else:
		protein_temp_anno_dir = temp_fasAnno + '/' + omaId
	species_blast_dir = hamstr_blast_dir + '/' + hamstrId + '/' + hamstrId
	#print species_blast_dir
	if not os.path.exists(hamstr_blast_dir + '/' + hamstrId):
		print 'ERROR: Blast database not found for species ', hamstrId
		sys.exit()
	# Get the ortholog (or seed protein) identifier present in the annotations directory
	temporary_file_1 = "temp_query_fas_%s.fa" %omaId
	temporary_file_2 = "temp_out_fas_%s.txt" %omaId
	ftemp = open(temporary_file_1, "w")
	ftemp.write('>' + omaId + '\n' + sequence)
	ftemp.close()
	blast_command = '%s -query %s -db %s -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out %s' %(blastp, temporary_file_1, species_blast_dir, temporary_file_2)
	os.system(blast_command)
	if os.path.exists(temporary_file_2) and len(open(temporary_file_2).read().split('\n')) > 1:
		hit_id = open(temporary_file_2).read().split('\n')[0].split('\t')[1]
		#print hit_id
		if omaId.split('_')[1] == species_id:
			fasCacheData.append(omaId.split('_')[1] + '\t' + hamstrId + '\t' + hit_id + '\n')
			logDict[omaId.split('_')[1]] = omaId.split('_')[1] + ':' + hit_id
		else:
			fasCacheData.append(omaId + '\t' + hamstrId + '\t' + hit_id + '\n')
			logDict[omaId] = omaId + ':' + hit_id
	if delTemp:
		os.system('rm -rf %s' %temporary_file_1)
		os.system('rm -rf %s' %temporary_file_2)

	# Extract FAS annotations for the protein sequence
	if not os.path.exists(protein_temp_anno_dir):
		anno_command = 'perl %s -path=%s -name=%s -extract=%s' %(annotation_script, species_anno_dir, hit_id, protein_temp_anno_dir)
		print 'Annotation command: ', anno_command				
		os.system(anno_command)
	return oma_ids, mapDict, fasCacheData, logDict 

### Actual FAS calculations - GreedyFAS is run from this module ###
def actual_FAS_exec(elements):
	fasFileData = []
	domainArchiData = []

	fas_score = "NA"
	if not elements == species_id:
		fas_command = 'python %s -s %s -p %s -r %s -o 1 -j temp_fas_score_%s.xml' %(greedy_fas_script, temp_fasAnno + '/' + species_id, temp_fasAnno + '/' + elements, speciesAnnoDir + '/' + mapDict[species_id], elements)				
		print 'FAS command: ', fas_command
		os.system(fas_command)

		# Parsing FAS output file
		print 'Parsing FACT score..'
		if os.path.exists('temp_fas_score_%s.xml' %elements):
			element = ET.parse('temp_fas_score_%s.xml' %elements).getroot()
			#print element
			for atype in element.findall('set_protein'):
				fas_score = atype.get('score')

		fasFileData.append(species_id + '\t' + logDict[species_id] + '\t' + elements + '\t' + logDict[elements] + '\t' + fas_score + '\n')
		
		print 'Extracting feature architecture for PhyloProfile..'
		domainArchiData = parseFAS('temp_fas_score_%s.xml' %elements, prot_id, species_id, elements)			

		#sys.exit()
		if delTemp:
			os.system('rm -rf temp_fas_score_%s.xml' %elements)

	#print fasFileData
	#print domainArchiData
	return fasFileData, domainArchiData

# FAS Calculation
def calculateFAS(working_dir, hamstr, spAnnoDir, orth_file, fas_file, domain_archi_file, map_file, blastp, protein_id, species_id, delTemp, hamstr_env):

	global prot_id
	prot_id = protein_id
	global speciesAnnoDir
	speciesAnnoDir = spAnnoDir

	# Using default paths for FAS scripts
	fas_dir = hamstr + '/bin/fas'

	global annotation_script
	annotation_script = fas_dir + '/annotation.pl'
	global greedy_fas_script
	greedy_fas_script = fas_dir + '/greedyFAS.py'

	### CHANGE HERE ###
	# Using default paths for HaMStR BLAST directory (where species gene sets are present)
	global hamstr_blast_dir
	if not hamstr_env == "":
		hamstr_blast_dir = hamstr + '/blast_dir_' + hamstr_env
	else:
		hamstr_blast_dir = hamstr + '/blast_dir'

	# Temporary directory to save FAS annotations for single sequences (seed and ortholog)
	global temp_fasAnno
	temp_fasAnno = working_dir + '/temp_fasAnnotations'

	if not os.path.exists(temp_fasAnno):
		os.mkdir(temp_fasAnno)

	# Read ortholog file into an array - Each element has 2 lines
	orthFileData = []

	orthFileRead = open(orth_file, "r")

	while True:
		line1 = orthFileRead.readline()
		line2 = orthFileRead.readline()
		if not line2:
			break
		orthFileData.append("\n".join([line1.split()[0], line2.split()[0]]))

	#print orthFileData
	#sys.exit()
	
	try:
		pool = Pool(processes=nr_processors)
		results = pool.map(retrieve_FAS_annotations, orthFileData)
	except KeyboardInterrupt as e:
		sys.exit(e)
	except:
		print("ERROR: Multiprocessing step <-> FAS annotations.")
		pass

	# Collecting data after multiprocessing run
	global mapDict, logDict
	oma_ids = []
	mapDict = {}
	fasCacheData = []
	logDict = {}

	for res in results:
		oma_ids.append(res[0][0])
		for key, value in res[1].iteritems():
			mapDict[key] = value
		fasCacheData.append(res[2][0])
		for key, value in res[3].iteritems():
			logDict[key] = value

	#print fasCacheData

	fasCacheFile = open("%s_cacheData.txt" %prot_id, "w")
	fasCacheFile.write("### Orthologs identifiers ###\n")
	for data in fasCacheData:
		fasCacheFile.write(data)
	fasCacheFile.close()
	

	# Calculate FAS scores between seed protein versus orthologs
	# Write the output in fas_file

	try:
		pool = Pool(processes=nr_processors)
		results = pool.map(actual_FAS_exec, oma_ids)
	except KeyboardInterrupt as e:
		sys.exit(e)
	except:
		print("ERROR: Multiprocessing step <-> FAS calculations.")
		pass

	fasFileData = []
	domainArchiData = []

	for res in results:
		if not res[0] == []:
			fasFileData.append(res[0][0])
		if not res[1] == []:
			for resListElements in res[1]:
				domainArchiData.append(resListElements)

	fasFile = open(fas_file, 'w')
	for data in fasFileData:
		fasFile.write(data)
	fasFile.close()

	domainFile = open(domain_archi_file, 'w')
	for data in domainArchiData:
		domainFile.write(data)
	domainFile.close()	

	if delTemp:
		os.system('rm -rf %s' %temp_fasAnno)

# HaMStROneSeq run
def run_hamstrOneSeq(hamstr, orth_file, map_file, prot_id, formatdb, blastp, proteome, delTemp, hamstr_env, include_paralogs):

	# Setting default paths for HaMStR-OneSeq
	hamstrOneSeq = hamstr + '/bin/oneSeq.pl'
	try:
		coreOrthDir = hamstr + '/core_orthologs'
		if not hamstr_env == "":
			taxaPath = hamstr + '/genome_dir_' + hamstr_env
		else:
			taxaPath = hamstr + '/genome_dir'	

		for line in open(orth_file):
			if line[0] == '>':
				omaId = line.split('\n')[0][1:6]
			else:
				protSeq = line.split('\n')[0]
		
		for line in open(map_file):
			if omaId in line.split()[-1]:
				hamstrId = line.split('\t')[0]
				break
		print 'hamstr id: ', hamstrId
		for taxas in glob.glob(taxaPath + '/*'):
			if hamstrId in taxas:
				refSpec = taxas.split('/')[-1]
				break
		
		with open(taxaPath + '/' + refSpec + '/' + refSpec + '.fa') as f:
			print 'Searching for the seqId..'
			flag = True
			for line in f:
				if f.next().split('\n')[0].replace('*', '').replace('X', '') == protSeq.replace('*', '').replace('X', ''):
					seqId = line.split('>')[1].split('\n')[0]
					flag = False
					break
		if flag:
			seqId = "NA"
			print 'No matching sequence found! Running BLAST search now..'
			#print 'Current working directory..'
			currentWorkDir = os.getcwd()
			#print 'Create a temporary directory..'
			tempDir = 'temp_blast_' + prot_id
			if not os.path.exists(tempDir):
				os.mkdir(tempDir)
			#print 'Change to the temporary directory..'
			os.chdir(tempDir)
			#print 'Create a temporary file with the input sequence..'
			#print 'Copy the reference proteome file into temporary directory..'
			os.system('cp -avr %s .' %(taxaPath + '/' + refSpec + '/' + refSpec + '.fa'))
			com = '%s -i %s' %(formatdb, refSpec +'.fa')
			#print 'Create blast database for the OMA sequences: ', com
			os.system(com)
			#print 'Perform BLAST search and pick up the top hit as query input ID..'
			os.system('%s -query %s -db %s -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp, orth_file, refSpec + '.fa'))
			#os.system('%s -query %s -db %s -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp, orth_file, proteome))
			if os.path.exists('temp.txt') and len(open('temp.txt').read().split('\n')) > 1:
				seqId = open('temp.txt').read().split('\n')[0].split('\t')[1]
			os.system('%s -query %s -db %s -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp, orth_file, proteome))
			if os.path.exists('temp.txt') and len(open('temp.txt').read().split('\n')) > 1:
				omaId = open('temp.txt').read().split('\n')[0].split('\t')[1]
			#print 'Writing OMA id into file'
			oma = open(omaIdFile, 'w')
			oma.write(omaId)
			oma.close()
			#print 'Change directory to original one..'
			os.chdir(currentWorkDir)
			#print 'Remove the temporary directory..'
			if delTemp:
				os.system('rm -rf %s' %tempDir)

		print 'SeqId: %s  ...  OmaId: %s' %(seqId, omaId)
		try:
			if not seqId == "NA" and delTemp:
				if os.path.exists('inputTaxaSet_oneSeq.txt'):
					if not hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -cleanup -fasoff -local -strict -rep -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif not hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -cleanup -fasoff -local -strict -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -cleanup -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
					elif hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -cleanup -fasoff -local -strict -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				else:
					if not hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -fasoff -local -strict -rep -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif not hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -fasoff -local -strict -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
					elif hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -fasoff -local -strict -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				
				print '##### Running hamstrOneSeq command: ', command
				os.system(command)

			elif not seqId == "NA" and not delTemp:
				if os.path.exists('inputTaxaSet_oneSeq.txt'):
					if not hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -fasoff -local -strict -rep -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif not hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -fasoff -local -strict -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
					elif hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -cleanup -coreTaxa="inputTaxaSet_oneSeq.txt" -fasoff -local -strict -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				else:
					if not hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -fasoff -local -strict -rep -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif not hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -fasoff -local -strict -seqName=%s -addenv=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id, hamstr_env)
					elif hamstr_env == "" and not include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -fasoff -local -strict -rep -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
					elif hamstr_env == "" and include_paralogs:
						command = 'perl %s -sequence_file=%s -seqid=%s -refSpec=%s -coreOrth=5 -minDist=species -maxDist=superkingdom -checkCoorthologsRef -fasoff -local -strict -seqName=%s' %(hamstrOneSeq, orth_file.split('/')[-1], seqId, refSpec, prot_id)
				
				print '##### Running hamstrOneSeq command: ', command
				os.system(command)

			else:
				print 'WARNING: No sequence id found to run HaMStR-OneSeq!!!'
		except:
			pass
			print 'HaMStROneSeq did not run properly!!!'
	except IOError:
		print 'WARNING: hamstrOneSeq did not run properly!!!'


	# Read the output files generated by OneSeq run (.fa and .fa.extended)
	# Write the results back to the original orthologs file
	print 'Adding HaMStR-OneSeq results to the orthologs file..'
	try:
		orthFileDefault = open(orth_file).read().split('\n')
		speciesList = {}
		speciesList[orthFileDefault[0][1:6]] = []
		speciesList[orthFileDefault[0][1:6]].append(orthFileDefault[1])
		fnew = open(orth_file, 'w')
		fnew.write(orthFileDefault[0] + '\n' + orthFileDefault[1] + '\n')

		#faFile = coreOrthDir + prot_id + '/%s.fa' %prot_id
		extendedFile = work_dir + '/%s.extended.fa' %prot_id
		'''if os.path.exists(faFile):
			omaId = "NA"
			with open(faFile) as f:
				for line in f:
					if line[0] == '>':
						hamstrId = line.split('|')[1].split('_')[0] + '_' + line.split('|')[1].split('_')[1]
						for m in open(map_file):
							if m.split('\t')[0] == hamstrId:
								omaId = m.split('\t')[3].split('\n')[0]
								break
						if not omaId == "NA" and omaId not in speciesList:
							speciesList.append(omaId)
							#fnew.write('>' + omaId + line.split('|')[2].split('\n')[0] + '\n' + f.next().replace('*', ''))
							fnew.write('>' + omaId + '\n' + f.next().replace('*', ''))
			#fnew.write('\n')'''
		if os.path.exists(extendedFile):
			#print 'Found .extended file..'
			omaId = "NA"
			with open(extendedFile) as f:
				for line in f:
					if '>' in line:
						hamstrId = line.split('|')[1]
						seqIdentifier = line.split('|')[2]
						seqHamstrOut = f.next().split()[0].replace('*', '').replace('X', '')
						for m in open(map_file):
							if hamstrId in m:
								omaId = m.split()[-1]
								break
						if not omaId == "NA":
							if omaId in speciesList.keys() and not seqHamstrOut in speciesList[omaId]:
								speciesList[omaId].append(seqHamstrOut)
								fnew.write('>' + str(len(speciesList[omaId]) + 1) + '_' + omaId + '\n' + seqHamstrOut + '\n')
							elif not omaId in speciesList.keys():
								speciesList[omaId] = []
								speciesList[omaId].append(seqHamstrOut)
								fnew.write('>1_' + omaId + '\n' + seqHamstrOut + '\n')
		fnew.close()
	except IOError:
		print 'ERROR: While writing HaMStr-OneSeq results!!!'

	if delTemp:
		if os.path.exists(coreOrthDir + '/' + prot_id):
			os.system('rm -Rf %s' %(coreOrthDir + '/' + prot_id))
		if os.path.exists("inputTaxaSet_oneSeq.txt"):
			os.system('rm inputTaxaSet_oneSeq.txt')


# Prepares input configuration file for REvolver
def prepareXML(xml_file, pfamDB, hmmfetch, aaMatrix, indel, p, sf, simTree, prot_id, hmm_file, output_dir):
	fnew = open(xml_file, 'w')
	fnew.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
	fnew.write('<configdata  xsi:schemaLocation="http://www.cibiv.at/Revolver ./input_schema.xsd" xmlns="http://www.cibiv.at/Revolver" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" >\n')
	fnew.write('\t<config>\n')
	fnew.write('\t\t<hmmdb path="%s"/>\n' %os.path.abspath(pfamDB))
	fnew.write('\t\t<hmmfetch location="%s"/>\n' %hmmfetch)
	fnew.write('\t</config>\n')
	fnew.write('\t<model>\n')
	fnew.write('\t\t<substitution name="%s"/>\n' %aaMatrix)
	fnew.write('\t\t<indel>\n')
	fnew.write('\t\t\t<insertion rate="%s">\n' %str(indel))
	fnew.write('\t\t\t\t<length distribution="geometric" p="%s"/>\n' %str(p))
	fnew.write('\t\t\t</insertion>\n')
	fnew.write('\t\t\t<deletion rate="%s">\n' %str(indel))
	fnew.write('\t\t\t\t<length distribution="geometric" p="%s"/>\n' %str(p))
	fnew.write('\t\t\t</deletion>\n')
	fnew.write('\t\t</indel>\n')
	fnew.write('\t</model>\n')
	fnew.write('\t<tree scalingFactor="%s" path="%s"  />\n' %(str(sf), os.path.abspath(simTree)))
	fnew.write('\t<root>\n')
	fnew.write('\t\t<inputSequence>\n')
	fnew.write('\t\t\t<fasta file="seq_%s.fa"/>\n' %prot_id)
	fnew.write('\t\t\t<hmmer file="%s"/>\n' %os.path.abspath(hmm_file))
	fnew.write('\t\t</inputSequence>\n')
	fnew.write('\t</root>\n')
	fnew.write('\t<output>\n')
	fnew.write('\t\t<dir path="%s" separateFastaFiles="false" trueAlignment="false" include="leaf"/>\n' %os.path.abspath(output_dir))
	fnew.write('\t</output>\n')
	fnew.write('</configdata>')
	
	fnew.close()

# Runs hmmscan and prepares the domain constraint file for REvolver
def hmmscan(hmmscan, orth_file, pfamDB, hmm_file, prot_id, species_id):
	ftemp = open('seq_%s.fa' %prot_id, 'w')
	with open(orth_file) as f:
		for line in f:
			if line[:3] == '>1_' and line.split()[0].split('_')[1] == species_id:
				ftemp.write(line.replace('_1', '') + f.next())
				break
	ftemp.close()
	os.system('%s --notextw -E 0.01 %s seq_%s.fa > %s' %(hmmscan, pfamDB, prot_id, hmm_file))
	#os.remove('tempFile.fa')
	

# Calculates indels rates
def calculateIndels(tree_file, trans, alnLength, iqtree24, def_indel, def_indel_dist):
	indel = float(def_indel)
	p = float(def_indel_dist)

	print '##### Calculating indels #####'
	try:
		trees = dendropy.TreeList.get_from_path(tree_file, "newick")
		tree_lengths = [tree.length() for tree in trees]
	except:
		pass
	
	result = ''
	try:
		command = "%s -s %s %s -tina -st MULTI" %(iqtree24, os.path.abspath(trans), os.path.abspath(tree_file))
		print 'IQ-Tree24 command: ', command
		result = subprocess.check_output(command, shell=True)
	except:
		print 'WARNING: IQTree-24 did not run properly!!!'
		pass
	
	for line in result.split('\n'):
		if line.split(':')[0] == 'mean length':
			if float(line.split(':')[1].replace(' ', '')) > 0:
				p = 1 / float(line.split(':')[1].replace(' ', ''))
				if p >= 1:
					p = 0.99
				elif p < 0.02:
					p = 0.02
			else:
				p = float(prot_config.default_indel_distribution)
		elif line.split(':')[0] == 'Parsimony score is':
			indel = (float(line.split(':')[1].replace(' ', '')) / (alnLength * tree_lengths[0])) / 2
	print 'Indel: ', indel

	fnew = open(indel_file, 'w')
	fnew.write(str(indel) + '\n' + str(p))
	fnew.close()
	

# Perform MSA of the ortholog sequences
# Convert the .aln format to .phy format
def performMSA(msa, clustalw):
	if cache:
		if os.path.exists(phy_file):
			print 'Re-using previously compiled orthologs alignment file'
		else:	
			try:
				os.system('%s %s > %s' %(msa, orth_file, aln_file))
				os.system('%s -convert -output=PHYLIP -infile=%s -outfile=%s' %(clustalw, aln_file, phy_file))
			except:
				pass
				print 'WARNING: MSA didn\'t work. Less than 2 sequences found for alignment!!!'
	else:
		try:
			os.system('%s %s > %s' %(msa, orth_file, aln_file))
			os.system('%s -convert -output=PHYLIP -infile=%s -outfile=%s' %(clustalw, aln_file, phy_file))
		except:
			pass
			print 'WARNING: MSA didn\'t work. Less than 2 sequences found for alignment!!!'

# Read OMA sequences file and parse OMA orthologs sequences
def findOmaSequences(prot_id, omaSeqs, species_id, mapFile):
	try:
		mapIds = []
		for line in open(mapFile):
			mapIds.append(line.split()[-1])

		print '##### Searching OMA orthologs sequences for %s #####' %prot_id
		fnew = open(orth_file, 'w')
		ids = open(id_file).read().split('\n')
		with open(omaSeqs) as f:
			for line in f:
				if line[0] == '>' and line.split('\n')[0][1:] in ids and line[1:6] in mapIds:
					#if line.split('\n')[0][1:6] == species_id:
						#fnew.write('>' + species_id + '\n' + f.next().replace('*', ''))
					fnew.write('>' + line[1:6] + '\n' + f.next().replace('*', '').replace('X', ''))
					#else:
					#	fnew.write(line + f.next().replace('*', ''))
		fnew.close()
	except IOError:
		sys.exit('ERROR: Cannot find OMA orthologs sequences. OMA sequence file does not exist!')
				
	
# Read OMA orthologs groups file and parses the ortholog list for input OMA id
def findOmaOrthologs(prot_id, querySeq, omaGroup, omaSeqs, proteome_file, formatdb, blastp, delTemp, species_id):
	try:
		if not querySeq == 'None' and speciesFoundInOmaFlag:
			run = 1
			sequenceFoundFlag = False
			with open(omaSeqs) as f:
				for line in f:
					if '>' in line and species_id in line:
						if f.next().split()[0].replace('*', "") == querySeq.split()[0].replace('*', ''):
							prot_id_temp = line.split('\n')[0][1:]
							oma = open(omaIdFile, 'w')
							oma.write(prot_id_temp)
							oma.close()
							#print 'Found OMA id:', inputId
							sequenceFoundFlag = True
							#fnew.write(inputId + '\n')
							#fnew.close()
							break

			if not sequenceFoundFlag:
				print '##### WARNING: The input sequence is not present in OMA sequences. Performing a BLAST search to get the closest best hit. #####'
				currentWorkDir = os.getcwd()
				tempDir = 'temp_blast_' + prot_id
				if not os.path.exists(tempDir):
					os.mkdir(tempDir)
				os.chdir(tempDir)
				os.system('cp -avr %s proteome.fa' %(proteome_file))
				com = '%s -i proteome.fa' %(formatdb)
				os.system(com)
				tempFile = open('temp_inputSeq.fa', 'w')
				tempFile.write('>' + prot_id + '\n' + querySeq)
				tempFile.close()
				os.system('%s -query temp_inputSeq.fa -db proteome.fa -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp))
				if os.path.exists('temp.txt') and len(open('temp.txt').read().split('\n')) > 1:
					prot_id_temp = open('temp.txt').read().split('\n')[0].split('\t')[1]
					oma = open(omaIdFile, 'w')
					oma.write(prot_id_temp)
					oma.close()
					sequenceFoundFlag = True
				os.chdir(currentWorkDir)
				if delTemp:
					os.system('rm -rf %s' %tempDir)

			if sequenceFoundFlag:
				run = 2
				print '##### Searching for OMA ortholog group for %s #####' %prot_id
				fnew = open(id_file, 'w')
				written = False
				for line in open(omaGroup):
					if not '#' in line and prot_id_temp in line:
						for ids in line.split('\n')[0].split('\t')[2:]:
							fnew.write(ids + '\n')
						fnew.close()
						written = True
						break
				if not written:
					run = 1
					fnew.write(prot_id_temp)
					fnew.close()
		elif not querySeq == 'None' and not speciesFoundInOmaFlag:
			run = 1
			print '##### The input sequence is not present in OMA database. Saving hit id from local blast directory. #####'
			currentWorkDir = os.getcwd()
			tempDir = 'temp_blast_' + prot_id
			if not os.path.exists(tempDir):
				os.mkdir(tempDir)
			os.chdir(tempDir)
			os.system('cp -avr %s proteome.fa' %(proteome_file))
			com = '%s -i proteome.fa' %(formatdb)
			os.system(com)
			tempFile = open('temp_inputSeq.fa', 'w')
			tempFile.write('>' + prot_id + '\n' + querySeq)
			tempFile.close()
			os.system('%s -query temp_inputSeq.fa -db proteome.fa -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp.txt' %(blastp))
			if os.path.exists('temp.txt') and len(open('temp.txt').read().split('\n')) > 1:
				prot_id_temp = open('temp.txt').read().split('\n')[0].split('\t')[1]
				oma = open(omaIdFile, 'w')
				oma.write(prot_id_temp)
				oma.close()
			os.chdir(currentWorkDir)
			if delTemp:
				os.system('rm -rf %s' %tempDir)
		else:
			print 'OMA ids given..'
			run = 2
			print '##### Searching for OMA ortholog group for already given OMA id %s #####' %prot_id
			oma = open(omaIdFile, 'w')
			oma.write(prot_id)
			oma.close()
			fnew = open(id_file, 'w')
			written = False
			for line in open(omaGroup):
				if not '#' in line and prot_id in line:
					for ids in line.split('\n')[0].split('\t')[2:]:
						fnew.write(ids + '\n')
					fnew.close()
					written = True
					break
			if not written:
				fnew.write(prot_id)
				fnew.close()
			
	except IOError:
		sys.exit('ERROR: Cannot find OMA orthologs id. Check OMA files given as input!')

	return run

# Read in OMA sequences file and create a new proteome file for the species id
# Performs formatdb on the proteome to be later used in reciprocal BLAST search
def parseProteome(species_id, omaSeqs, makeblastdb, proteome_file, crossRefFile, hamstrDir, hamstrEnv, search_oma_database):
	try:
		# Check whether the parsed proteome in already present in Cache directory
		if cache and os.path.exists(cache_dir + '/proteome_' + species_id):
			print '##### Gene set for species %s found in Cache directory. Reusing it. #####' %species_id
			os.system('cp -avr %s %s' %(cache_dir + '/proteome_' + species_id, proteome_file))
		elif search_oma_database:
			print '##### Parsing gene set for species %s from OMA database #####' %species_id
			fnew = open(proteome_file, 'w')
			
			speciesFoundInOmaFlag = False
			with open(omaSeqs) as f:
				for line in f:
					if line[:6] == '>' + species_id:
						speciesFoundInOmaFlag = True
						fnew.write(line + f.next())
			fnew.close()

			if not speciesFoundInOmaFlag:
				sys.exit('ERROR: Species %s not found in OMA database. Please make sure if the species exists in OMA database. If not, please turn off "search_oma_database" flag in program configuration file.' %species_id)

		else:
			print('##### Parsing gene set for species %s in local (HaMStR) blast directory. #####' %species_id)
			for line in open(crossRefFile):
				if species_id == line.split()[-1]:
					hamstr_id = line.split()[0]
					break

			if not hamstrEnv == "":
				blast_dir = hamstrDir + '/blast_dir_' + hamstrEnv
			else:
				blast_dir = hamstrDir + '/blast_dir'

			species_blast_file = blast_dir + '/' + hamstr_id + '/' + hamstr_id + '.fa'
			
			if os.path.exists(species_blast_file):
				os.system('cp -avr %s %s' %(species_blast_file, proteome_file))
			else:
				sys.exit('ERROR: Reference species %s not found in local (HaMStR) BLAST directory!!!' %species_id)
			
			os.system('cp -avr %s %s' %(proteome_file, cache_dir + '/proteome_' + species_id))
		print '##### Making BLAST db of the gene set to be used by the blast search #####'
		os.system('%s -in %s -input_type fasta -dbtype prot' %(makeblastdb, proteome_file))
		#proteome_file = os.path.abspath(proteome_file)

	except IOError:
		sys.exit('ERROR: Cannot create gene set for species %s. Please check the path of species information contatining file!' %species_id)
