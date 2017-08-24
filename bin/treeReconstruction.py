import os, sys
import subprocess
import time
import maxLikDistMatrix

# Module to reconstruct tree using RAxML and do scaling factor calculation

# Calculate the median of a set
def median(lst):
    even = (0 if len(lst) % 2 else 1) + 1
    half = (len(lst) - 1) / 2
    return sum(sorted(lst)[half:half + even]) / float(even)

# Rename the ortholog group sequences into oma ids
def rename_orth_file():
	fnew = open('temp_orth_%s.fa' %protein_id, 'w')
	for i in range(0, len(orth_file) - 1, 2):
		fnew.write(orth_file[i][:6] + '\n' + orth_file[i+1] + '\n')
	fnew.close()

# Perform MSA of the new renamed ortholog sequences file (MAFFT linsi)
# Perform degapping using in-house script by Ingo (degapper.pl)
# Convert the MSA format from .aln to .phy (clustalw2)
def msa_convert():
	#print 'MSA of the renamed ortholog sequences file:'
	if os.path.exists(aln_file):
		os.system('cp -avr %s temp_orth_%s.aln' %(aln_file, protein_id))
	else:
		os.system('%s temp_orth_%s.fa > temp_orth_%s.aln' %(linsi, protein_id, protein_id))

	if os.path.exists(phy_file):
		os.system('cp -avr %s temp_orth_%s.phy' %(phy_file, protein_id))
	else:
		os.system('%s -convert -infile=temp_orth_%s.aln -outfile=temp_orth_%s.phy -output=PHYLIP' %(clustalw, protein_id, protein_id))

	##### ACTIVATE THIS BLOCK IF DEGAPPING HAS TO BE IMPLEMENTED #####
	'''#print 'File format conversion (.aln -> .fa):'
	s1 = open('temp_orth_%s.aln' %protein_id).read().split('\n')
	fnew = open('temp_orth_aln_%s.fa' %protein_id, 'w')
	fnew.write(s1[0] + '\n')
	for i in range(1, len(s1) - 1):
		if s1[i][0] == '>':
			fnew.write('\n' + s1[i] + '\n')
		else:
			fnew.write(s1[i])
	fnew.close()
	#print 'Apply degapping algorithm (degapper.pl):'
	os.system('perl %s -limit=0.50 -in=temp_orth_aln_%s.fa -out=temp_orth_aln_degap_%s.fa' %(degap, protein_id, protein_id))
	
	#print 'File format conversion (.fa -> .phy):'
	os.system('%s -convert -infile=temp_orth_aln_degap_%s.fa -outfile=temp_orth_%s.phy -output=PHYLIP' %(clustalw, protein_id, protein_id))
	#print 'complete..'
	'''

# Run RAxML for tree reconstruction
def run_raxml():
	if makeTree:
		if reuse_cache and os.path.exists('RAxML_bestTree.' + protein_id):
			print('Orthologs RAxML tree found. Reusing it.')
		else:
			os.system('rm -rf RAxML_*')
			os.system('%s-PTHREADS -T %s -s temp_orth_%s.phy -m PROTGAMMA%s -p 12345 -n %s' %(raxml, nr_proc, protein_id, aaMatrix, protein_id))
	#print 'complete..'

# Rename the reconstructed tree file
def rename_raxml():
	try:
		s1 = open('RAxML_bestTree.' + protein_id).read()
		fnew = open('RAxML_bestTree.' + protein_id, 'w')
		for i in range(0, len(orth_file) - 1, 2):
			nick_name = orth_file[i][1:6]
			full_name = orth_file[i][1:]
			s1 = s1.replace(nick_name, full_name)
		fnew.write(s1)
		fnew.close()
	except IOError:
		sys.exit('RAxML tree could not be found!')

# Remove all the temp files generated
def rm_temp():
	
	os.system('rm temp_orth_%s.fa' %protein_id)
	os.system('rm temp_orth_%s.aln' %protein_id)
	os.system('rm temp_orth_%s.phy' %protein_id)
	os.system('rm temp_orth_%s.phy*' %protein_id)
	if os.path.exists('temp_orth_%s.phy.reduced' %protein_id):
		os.system('rm temp_orth_%s.phy.reduced' %protein_id)
	os.system('rm temp_parameters_%s.txt' %protein_id)
	os.system('rm maxLikDist_%s.txt' %protein_id)
	os.system('rm RAxML_log.*')
	os.system('rm RAxML_parsimonyTree*')
	os.system('rm RAxML_result*')
	os.system('rm RAxML_info*')

# Perform likelihood mapping on the aligned sequences (for 4 or more sequences)
def likelihoodMapping():
	
	fnew = open('temp_parameters_%s.txt' %protein_id, 'w')
	# Select the parameter file for puzzle i.e. fill in the alignmentFile used !!! *** CHANGE THE FILE HERE ***
	s1 = open(params_puzzle).read()
	fnew.write(s1.replace('alignmentFile', 'temp_orth_%s.phy' %protein_id))
	fnew.close()

	#print 'Running tree puzzle..'
	os.system('%s < temp_parameters_%s.txt' %(puzzle,protein_id))

	#print 'Puzzle run complete..'

# Perform likelihood mapping on the aligned sequences (for 3 or more sequences)
def distanceMapping():

	fnew = open('temp_parameters_%s.txt' %protein_id, 'w')
	# Select the parameter file for puzzle i.e. fill in the alignmentFile used !!! *** CHANGE THE FILE HERE ***
	s1 = open(params_puzzle).read()
	fnew.write(s1.replace('alignmentFile', 'temp_orth_%s.phy' %protein_id).replace('b', 'k\nk'))
	fnew.close()

	#print 'Running tree puzzle..'
	os.system('%s < temp_parameters_%s.txt' %(puzzle, protein_id))
	
	#print 'Puzzle run complete..'

# Calculate the scaling factor based on maximum likelihood distances
def scalingFactorMax():
	scales = []
	# Generate maximum likelihood distance file for orthologs
	### NOTE THE FILE USED HERE!!!!*******************************************!!!!
	outfile = open('temp_orth_%s.phy.dist' %protein_id).read().split('\n')
	maxLikDistMatrix.main(outfile, protein_id)
	try:
		orthMaxFile = open('maxLikDist_%s.txt' %protein_id).read().split('\n')
		speciesMaxFile = open(species_maxLikMatrix).read().split('\n')
		hamstrFile = open(map_file).read().split('\n')
		for i in range(len(orthMaxFile) - 1):
			line = orthMaxFile[i].split('\t')[0]
			if '_' in line:
				species1 = line.split('_')[1]
			else:
				species1 = line[1:]

			###
			### This block is needed when max likelihood matrix do not have OMA identifiers ###
			###
			#print species1
			#for j in range(len(hamstrFile) - 1):				
			#	if species1 == hamstrFile[j].split('\t')[3]:
			#		hamstr1 = hamstrFile[j].split('\t')[0]
			#		break

			for k in range(i + 1, len(orthMaxFile) - 1):
				line = orthMaxFile[k].split('\t')[0]
				if '_' in line:
					species2 = line.split('_')[1]
				else:
					species2 = line[1:]
				
				#for j in range(len(hamstrFile) - 1):
				#	if species2 == hamstrFile[j].split('\t')[3]:
				#		hamstr2 = hamstrFile[j].split('\t')[0]
				#		break
				#print species1, species2
				maxDistOrth = float(orthMaxFile[i].split('\t')[k + 1])
				maxDistSpecies = 0 #Maximum likelihood distance from species max. likelihood matrix
				#print maxDistOrth
				flag1 = True
				flag2 = True
				for l in range(len(speciesMaxFile) - 1):
					if speciesMaxFile[l].split('\t')[0] == species1:
						rowIndex = l
						flag1 = False
					elif speciesMaxFile[l].split('\t')[0] == species2:
						columnIndex = l
						flag2 = False
				#if not flag1 and not flag2:
				#	maxDistSpecies = float(speciesMaxFile[rowIndex].split('\t')[columnIndex])
				#else:
				#	maxDistSpecies = 1.0
				mlPresent = True
				if flag1 or flag2:
					#print 'yes'
					# Checking for the likelihood score in cache directory
					if os.path.exists(cacheDir + '/' + species1 + '_' + species2 + '.lik'):
						maxDistSpecies = float(open(cacheDir + '/' + species1 + '_' + species2 + '.lik').read().split('\n')[0])
					elif os.path.exists(cacheDir + '/' + species2 + '_' + species1 + '.lik'):
						maxDistSpecies = float(open(cacheDir + '/' + species2 + '_' + species1 + '.lik').read().split('\n')[0])
					else:
						#print 'No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2)
						#maxDistSpecies = 1.00
						mlPresent = False
				else:
					if not speciesMaxFile[rowIndex].split('\t')[columnIndex] == "NA":
						maxDistSpecies = float(speciesMaxFile[rowIndex].split('\t')[columnIndex])
					else:
						# Checking for the likelihood score in cache directory
						if os.path.exists(cacheDir + '/' + species1 + '_' + species2 + '.lik'):
							maxDistSpecies = float(open(cacheDir + '/' + species1 + '_' + species2 + '.lik').read().split('\n')[0])
						elif os.path.exists(cacheDir + '/' + species2 + '_' + species1 + '.lik'):
							maxDistSpecies = float(open(cacheDir + '/' + species2 + '_' + species1 + '.lik').read().split('\n')[0])
						else:
							#print 'No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2)
							mlPresent = False
							#maxDistSpecies = 1.00

				#print maxDistSpecies
				if mlPresent and not maxDistSpecies == 0:
					scales.append(maxDistOrth / maxDistSpecies)
				elif mlPresent and maxDistSpecies == 0:
					#scales.append(1.00)
					pass
	except:
		print '### ERROR: Scaling factor calculation had an error ###'
		sys.exit('Maximum likelihood files are invalid!')

	if len(scales) >= 1:
		#scale_value = sum(scales) / len(scales)
		#return scale_value
		#print 'Scales: ', scales
		#print 'Mean of Scales: ', sum(scales) / float(len(scales))
		return median(scales)
	else:
		return sf							

# Main module for running tree reconstruction
def main(Raxml, Linsi, Clustalw, Orthologs, AaMatrix, Protein_id, Puzzle, Params_puzzle, Map_file, Species_maxLikMatrix, Scale_file, Tree_file, delTemp, defScale, cache_dir, ortholog_tree_reconstruction, nr_processors, use_cache):

	global raxml, linsi, clustalw, orth_file, aaMatrix, protein_id, puzzle, params_puzzle, map_file, species_maxLikMatrix, scaleFile, treeFile, phy_file, aln_file, sf, cacheDir, makeTree, nr_proc, reuse_cache
	cacheDir = cache_dir
	raxml = Raxml
	linsi = Linsi
	clustalw = Clustalw
	aaMatrix = AaMatrix
	protein_id = Protein_id
	puzzle = Puzzle
	params_puzzle = Params_puzzle
	map_file = Map_file
	species_maxLikMatrix = Species_maxLikMatrix
	orth_file = open(Orthologs).read().split('\n')
	scaleFile = Scale_file
	treeFile = Tree_file
	aln_file = 'ogSeqs_' + protein_id + '.aln'
	phy_file = 'ogSeqs_' + protein_id + '.phy'
	sf = float(defScale)
	makeTree = ortholog_tree_reconstruction
	nr_proc = nr_processors
	reuse_cache = use_cache


	#sf = 1.00
		
	if len(orth_file) > 7:
		try:
			rename_orth_file()
			msa_convert()
			run_raxml()
			if reuse_cache and os.path.exists(scaleFile):
				print('Pre-computed scaling factor found. Reusing it.')
			else:
				likelihoodMapping()
				sf = scalingFactorMax()
				print 'Scaling factor: ', sf
				fnew = open(scaleFile, 'w')
				fnew.write(str(sf))
				fnew.close()
			
			if delTemp:
				rm_temp()
		except:
			print '### ERROR: Some step in the tree reconstruction was invalid!! ###'
			pass
		
					
		
	elif len(orth_file) < 8 and len(orth_file) > 3:
		try:
			rename_orth_file()
			msa_convert()
			distanceMapping()
			sf = scalingFactorMax()
			
		except:
			pass
		print 'Scaling factor: ', sf
		fnew = open(scaleFile, 'w')
		if sf > 0:
			fnew.write(str(sf))
		else:
			fnew.write(str(defScale))
		fnew.close()
		if delTemp:
			rm_temp()
	
