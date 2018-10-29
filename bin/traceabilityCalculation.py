import os, sys
import configure
import subprocess
import dendropy
import time
from multiprocessing import Pool

### Module where actual REvolver / BLAST cycles run takes place
### After every REvolver run, detection probability is checked with reciprocal Blast search and stored in a hash table
### Using the final hash table, decay rate is calculated using decay script

#
# Version history:
#	1. Arpit
#	2. Dominik (Bugfixes, changelist below)
#		- Multiprocessing (utilized in FAS annotation and FAS score calculation) did not receive any exceptions (a bug in python). Now, the main thread handles the CRTL-C induced interruption exception itself by calling pool.terminate and pool.join.
#######################################

def main(p_id, config_file):
	global prot_id
	prot_id = p_id

	rootDir = os.getcwd()

	global prot_config
	prot_config = configure.setParams(config_file)
	cache = prot_config.reuse_cache
	nr_proc = prot_config.nr_processors

	work_dir = prot_config.path_work_dir + '/' + prot_id
	os.chdir(work_dir)

	if os.path.exists(work_dir + '/omaId.txt'):
		global blastHitId
		blastHitId = open(work_dir + '/omaId.txt').read().split('\n')[0]
	else:
		sys.exit('### ERROR: No reciprocal BLAST hit id found!!!! ###')

	global xml_file
	xml_file = 'revolver_config_' + prot_id + '.xml'

	global proteome_file
	proteome_file = 'proteome_' + prot_id

	trees = dendropy.TreeList.get_from_path(prot_config.simulation_tree, "newick")
	global taxonset
	taxonset = []
	for element in trees.taxon_namespace:
		taxonset.append(str(element).replace("'", ""))
	taxonset = taxonset[::-1]

	print '##### Running REvolver / BLAST cycles: #####'
	start_time = time.time()

	if cache and os.path.exists('decay_summary_%s.txt_parameter' %prot_id):
		pass
	else:
		try:
			pool = Pool(processes=nr_proc)
			results = pool.map(actual_traceability_calculation, range(prot_config.simulation_runs))
		except KeyboardInterrupt as e:
			pool.terminate()
			pool.join()
			print("Interrupting REvolver")
			sys.exit(e)
		except:
			print("ERROR: Multiprocessing step <-> Traceability Calculations.")
			pass

		print '#####\tTIME TAKEN: %s mins REvolver/BLAST#####' %((time.time() - start_time) / 60)

		ffull = open('full_decay_results_%s.txt' %prot_id, 'w')
		fsum = open('decay_summary_%s.txt' %prot_id, 'w')

		detection_probability = {}
		for res in results:
			for key, value in res.iteritems():
				if not key in detection_probability.keys():
					detection_probability[key] = []
					detection_probability[key].append(value)
				else:
					detection_probability[key].append(value)
		#print detection_probability

		for taxa in taxonset:
			ffull.write(taxa + ' ')
			count = 0
			for element in detection_probability[taxa]:
				ffull.write(str(element))
				count += int(element)
			ffull.write('\n')
			fsum.write(str(float(count) / float(prot_config.simulation_runs)) + '\n')
		ffull.close()
		fsum.close()

		print '##### Calculating decay parameters #####'
		decayParams(prot_config.R, prot_id, prot_config.decay_script)

	os.chdir(rootDir)

def actual_traceability_calculation(run):
	temp_revolver_config_file = "%s_revolver_config.xml" %run
	os.system('cp %s %s' %(xml_file,temp_revolver_config_file))
	temp_data = open(temp_revolver_config_file).read()
	ftemp = open(temp_revolver_config_file, 'w')
	ftemp.write(temp_data.replace('REvolver_output', 'REvolver_output_' + str(run)))
	ftemp.close()

	revolver_output_dir = "REvolver_output_%s" %run
	if not os.path.exists(revolver_output_dir):
		os.mkdir(revolver_output_dir)

	detection_probability = {}

	print 'Run: ', run
<<<<<<< HEAD
=======
	#command = 'java -Xmx2G -Xms2G -cp "%s" revolver %s' %(prot_config.REvolver, temp_revolver_config_file)
	#print 'REvolver calculations command: ', command
>>>>>>> 1b295e1fcfd8cd322e31178a2e9e93e0f4d95144

	success = False
	trials = 0
	while(not success and trials < 10):
		trials += 1
		try:
			run_revolver(prot_config.REvolver, temp_revolver_config_file)
			blastOutput = run_blast(prot_config.blastp, prot_id, proteome_file, revolver_output_dir)
			for taxa in taxonset:
				detection = 0
				linecount = 1
				for line in blastOutput.split('\n'):
					#print line
					if taxa == line.split('\t')[0]:
						if line.split('\t')[1] == blastHitId:
							detection = 1
							break
						linecount += 1
						if linecount > 5:
							break
				detection_probability[taxa] = detection
			success = True
		except KeyboardInterrupt:
			print('Keyboard interruption by user!')
			raise Exception
		#except:
		#	pass
	if trials >= 10:
		print('TOO MANY TRIALS FOR REVOLVER!!! Check REvolver configuration file.')

	if prot_config.delete_temp:
		os.system('rm -rf %s' %revolver_output_dir)
		os.system('rm -rf %s' %temp_revolver_config_file)

	return detection_probability

def decayParams(r, prot_id, decay_script):
<<<<<<< HEAD
	command = '%s --quiet --vanilla --file=%s --args decay_summary_%s.txt' %(r, decay_script, prot_id)
=======
	command = '%s --quiet --vanilla %s decay_summary_%s.txt' %(r, decay_script, prot_id)
>>>>>>> 1b295e1fcfd8cd322e31178a2e9e93e0f4d95144
	print '##### Decay parameter calculation command: ', command
	os.system(command)

def run_revolver(REvolver, xml_file):
	command = 'java -Xmx2G -Xms2G -cp "%s" revolver %s' %(REvolver, xml_file)
	#print '##### REvolver calculations command: ', command
	os.system(command)

def run_blast(blastp, prot_id, proteome, revolverOut):
	command = '%s -query %s/out.fa -db %s -outfmt "6 qseqid sseqid"' %(blastp, revolverOut, proteome)
	result = subprocess.check_output(command, shell=True)
	return result

