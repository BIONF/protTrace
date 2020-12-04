import os, sys
import glob

### Supporting script for maximum likelihood calculations ###

### Orchestrates the calculation of species distances between 
### the query species and all other species in the 
### speciesTreeMapping list. This function can be called as 
### a __main__ to distribute ProtTrace tasks separately.
def CalculateSpeciesDistances(query,config):

    # Move to the dedicated distance calculation root 
    # directory for all species
    rootDir = config.path_distance_work_dir
    os.chdir(rootDir) 
    
    speciesMapping = config.hamstr_oma_tree_map

    # Gather the set of all subject species in this project 
    # from the speciesTreeMapping file
    speciesSet = set()
    try:
        with open(speciesMapping, 'r') as sm:
            for line in sm:
                speciesSet.add(line.split()[-1]
    except IOError:
        sys.exit('ERROR: Could not open %s. Please check the path.' %crossRefFile)

    # Exclude all species with existing distances in speciesMaxLikelihood
    missingSpecies = CheckSpeciesMaxLikelihoodTable(query,speciesSet,config)

    # Exclude missing species with existing .lik files in the cache directory
    if len(missingSpecies) > 0:
        os.chdir(
    
# Checks the speciesMaxLikelihood file for existing distances
# Returns missing target species OMA IDs
def CheckSpeciesMaxLikelihoodTable(query,targets,config):

    # Check whether all query-subject pairs have already
    # computed distances in the speciesMaxLikelihood file
    speciesMaxLikelihood = config.species_MaxLikMatrix
    # We need a list first, because we first fetch the names
    # in the first line, then we fetch existing numbers
    # in the query species row, where both indices must
    # coalign
    computedSpeciesList = list()
    computedSpeciesSet = set()
    try:
        with open(speciesMaxLikelihood, 'r') as ml:
            # The first row contains all species names
            # Its first column is empty
            # We populate the species list with it
            # The file reader will be positioned at the next line anyways
            computedSpeciesList = ml.readline().rstrip().split("\t")[1:]
            # Now we look up our query species row for which distances are
            # actually computed / not N/A / not missing
            for row in ml:
                splitRow = row.split("\t")
                # Find the row where the first item corresponds to our query species
                if splitRow[0] == query:
                    # The first item was only necessary to recognize the query species
                    # This way, we save a +1 index operation for every column 
                    # (and the headache)
                    splitRow.Remove(0)
                    # Create a set of all species names whose index in this column 
                    # contains a valid digit
                    computedSpeciesSet = (computedSpeciesList[i] for i in range(len(computedSpeciesList)) if isdigit(splitRow[i]))
                    # We do not need the list anymore
                    computedSpeciesList = None
                    break
    except IOError:
        sys.exit('ERROR: Could not open %s. Please check the path.' %crossRefFile)
 
    # This is done to remove references and in case 
    # targets is a list
    targetSpeciesSet = set(targets)
    
    # Returns species in the speciesTreeMapping list without
    # digits in SpeciesMaxLikelihood
    return targetSpeciesSet.difference(computedSpeciesSet)

def CalculateProteinDistances(species1,species2,config):

    # Read the configuration of ProtTrace for paths
    omaSeqs = config.path_oma_seqs
    omaPairs = config.path_oma_pairs
    concatAlignment = config.concat_alignments_script
    omaProteomesDir = ".."
    linsi = config.msa
    clustalw = config.clustalw
    deleteTemp = config.delete_temp
    
    # Create a logging file
    if not os.path.exists('log.txt'):
        log = open('log.txt', 'w')
    else:
        log = open('log.txt', 'a')
    
    # Create the working directory for the processed species pair
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
    
    # Define the preprocessing step, where the set of orthologous 
    # protein alignments are gathered
    def preprocess():
        print('Preprocessing:\tParsing sequence pairs and aligning them...')
        flag = True
        c = 1
        for line in open(omaPairs):
            if species1 in line and species2 in line:
                flag = False
                print('Sequence Nr.:', c)
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
                    print('Searching in file: ', omaSeqs, species1)
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
                    print('Searching in file: ', omaSeqs, species2)
                    with open(omaSeqs) as omaseq:
                        for line2 in omaseq:
                            if '>' in line2 and species2Id in line2:
                                fnew.write(line2[:6] + '\n' + omaseq.next().replace('*', ''))
                                break
                fnew.close()
                c += 1
                        # ProtTrace rather uses MAFFT linsi than muscle
                        # To avoid adding another dependency, I uncommented the linsi command
                #os.system('muscle -quiet -in %s -out %s' %(filename, filename.replace('.fa', '.aln').replace(faDir, alnDir)))
                        os.system('%s %s > %s' %(linsi, filename, filename.replace('.fa', '.aln').replace(faDir, alnDir)))
    
        if flag:
            print('WARNING! ERROR! Species missing!')
            log.write('WARNING: Either %s or %s is missing in the OMA database' %(species1, species2))
    # Collect all pairwise protein alignments, concatenate them
    # and calculate the summarized pairwise species distance
    def postprocess():
        print('Preprocessing complete..\nConcating the alignments..')
        concatFile = species1+'_'+species2+'.aln'
        # The concatAlignment script requires perl
        os.system('perl %s -in=%s -out=%s' %(concatAlignment, alnDir, concatFile))
        print('Postprocessing concatenated alignment file..')
        postProcessConcatFile = concatFile.replace('.aln', '.dup.aln')
        fnew = open(postProcessConcatFile, 'w')
        with open(concatFile) as con:
            for line in con:
                if '>' in line:
                    seqId = line
                    sequence = con.next()
                    fnew.write(seqId + sequence + seqId.replace('>', '>dup_') + sequence)
        fnew.close()
        print('Converting alignment file to phylip format..')
        phyFile = postProcessConcatFile.replace('.aln', '.phy')
        os.system('%s -convert -output=phylip -infile=%s -outfile=%s' %(clustalw, postProcessConcatFile, phyFile))
    
        print('Performing likelihood mapping..')
        puzzleParams = 'temp_puzzleParams.txt'
        fnew = open(puzzleParams, 'w')
        fnew.write(phyFile + '\nb\ne\nm\nm\nm\nm\nm\nm\ny\n')
        fnew.close()
        os.system('puzzle < %s' %puzzleParams)
    
        print('Parsing outfile and saving likelihood distance to')
        result = species1+'_'+species2+'.lik'
        fnew = open(result, 'w')
        outdist = open(phyFile + '.dist').read().split('\n')
        fnew.write(outdist[3].split()[1])
        fnew.close()
    #if deleteTemp: 
    # Deletes temporary files 
    #    if os.path.exists(result) and not len(open(result).read().split('\n')) == 0:
    #        os.system('rm -rf %s' %faDir)
    #        os.system('rm -rf %s' %alnDir)
    #        os.system('rm %s*' %phyFile)
    #        os.system('rm -rf *.txt')
    #        os.system('rm -rf *.aln')

    # END OF POSTPROCESS FUNCTION DEFINITION
    
    # Performs the functions defined before
    os.chdir(workDir)
    preprocess()
    postprocess()
    os.chdir(rootDir)

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        print('ERROR:\tNo arguments entered for calculating species distances:\nUSAGE:\tdistanceCalculation.py -s <query species ID> -c <configFile> [-help]')
        sys.exit(2)
