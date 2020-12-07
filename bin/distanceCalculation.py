import os, sys
import glob

### Supporting script for maximum likelihood calculations ###

### Orchestrates the calculation of species distances between 
### the query species and all other species in the 
### speciesTreeMapping list. <config> is expected as an 
### instance of the ProtTrace configuration (generated with 
### configure.setParams)
def calculate_species_distances(config):

    # Read the query species from the ProtTrace configuration
    query = config.species

    # Move to the dedicated distance calculation root 
    # directory for all species
    root_dir = config.path_distance_work_dir
    os.chdir(root_dir)

    # Check and prepare the cache location
    cache_dir = config.path_cache
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
    
    species_mapping = config.hamstr_oma_tree_map

    # Gather the set of all subject species in this project 
    # from the speciesTreeMapping file
    species_set = set()
    try:
        with open(species_mapping, 'r') as sm:
            for line in sm:
                species_set.add(line.split()[-1])
    except IOError:
        sys.exit('ERROR: Could not open %s. Please check the path.' %crossRefFile)

    # Exclude all species with existing distances in species_max_likelihood
    missing_species = check_species_max_likelihood_table(query,species_set,config)

    # Exclude missing species with existing .lik files in the cache directory
    if len(missing_species) > 0:
        os.chdir(cache_dir)
        # We check the cache directory for {query species}_{missing species}.lik files
        missing_species = {species for species in missing_species if not os.path.isfile("{0}_{1}.lik".format(query,species))}
        os.chdir(root_dir)

    # If we are still missing species distances, we calculate them now
    # All distances are copied to the cache directory
    # This means they stay backed up when the not saving cache option would wipe them
    if len(missing_species) > 0:
        for species in missing_species:
             calculate_protein_distances(query,species,config,cache_dir,1)
             break
    
# Checks the species maximum likelihood file for existing distances
# Returns missing target species OMA IDs
def check_species_max_likelihood_table(query,targets,config):

    # Check whether all query-subject pairs have already
    # computed distances in the species_max_likelihood file
    species_max_likelihood = config.species_MaxLikMatrix
    # We need a list first, because we first fetch the names
    # in the first line, then we fetch existing numbers
    # in the query species row, where both indices must
    # coalign
    computed_species_list = list()
    computed_species_set = set()
    try:
        with open(species_max_likelihood, 'r') as ml:
            # The first row contains all species names
            # Its first column is empty
            # We populate the species list with it
            # The file reader will be positioned at the next line anyways
            computed_species_list = ml.readline().rstrip().split("\t")[1:]
            # Now we look up our query species row for which distances are
            # actually computed / not N/A / not missing
            for row in ml:
                split_row = row.rstrip().split("\t")
                # Find the row where the first item corresponds to our query species
                if split_row[0] == query:
                    # The first item was only necessary to recognize the query species
                    # This way, we save a +1 index operation for every column 
                    # (and the headache)
                    split_row = split_row[1:]
                    # Create a set of all species names whose index in this column 
                    # contains a valid digit
                    for i in range(len(computed_species_list)):
                        if split_row[i] is not None:
                            if split_row[i].isdigit():
                                computed_species_set.add(computed_species_list[i])
#                    computed_species_set = (computed_species_list[i] for i in range(len(computed_species_list)) if split_row[i] is not None and split_row[i].isdigit())
                    break
    except IOError:
        sys.exit('ERROR: Could not open %s. Please check the path.' %crossRefFile)
 
    # This is done to remove references and in case 
    # targets is a list
    target_species_set = set(targets)
    
    # Returns species in the speciesTreeMapping list without
    # digits in SpeciesMaxLikelihood
    return target_species_set.difference(computed_species_set)

### Calculates the pairwise species maximum likelihood distance 
### between species1 and species2. The distance is copied to the
### target_dir. The config is loaded from ProtTrace's configure.py
### if this script is executed directly.
def calculate_protein_distances(species1,species2,config,target_dir,add_filename=1):

    print("Calculating the protein distance between {0} and {1}.".format(species1,species2))

    # Read the configuration of ProtTrace for paths
    oma_seqs = config.path_oma_seqs
    oma_pairs = config.path_oma_pairs
    concatAlignment = config.concat_alignments_script
    oma_proteomes_dir = ".."
    linsi = config.msa
    #clustalw = config.clustalw
    nr_processors = config.nr_processors
    delete_temp = config.delete_temp
    
    # Create the working directory for the processed species pair
    work_dir = species1+'_'+species2
    
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    
    # The directory which contains the fasta files
    fa_dir = os.path.abspath(work_dir) + '/fa_dir'
    # The directory which contains the alignments
    aln_dir = os.path.abspath(work_dir) + '/aln_dir'
    
    if not os.path.exists(fa_dir):
        os.mkdir(fa_dir)
    
    if not os.path.exists(aln_dir):
        os.mkdir(aln_dir)

    os.chdir(work_dir)
    
    # Align the set of pairwise oprthologs between the query
    # and the target species
    print('Preprocessing:\tParsing sequence pairs and aligning them...')
    c = 1
    with open(oma_pairs,'r') as op:
        for line in op:
            if species1 in line and species2 in line:
                print('Sequence Nr.:', c)
                if species1 in line.split('\t')[0]:
                    species1Id = line.split('\t')[0]
                    species2Id = line.split('\n')[0].split('\t')[1]
                else:
                    species2Id = line.split('\t')[0]
                    species1Id = line.split('\n')[0].split('\t')[1]
    
                filename = fa_dir + '/seq_'+str(c)+'.fa'
                with open(filename, 'w') as fnew:
                    if os.path.exists(oma_proteomes_dir + '/proteome_' + species1):
                        with open(oma_proteomes_dir + '/proteome_' + species1) as prot1:
                            for prot_line in prot1:
                                if '>' in prot_line and species1Id in prot_line:
                                    fnew.write(prot_line[:6] + '\n' + prot1.next().replace('*', ''))
                                    break
    
                    else:
                        print('Searching in file: ', oma_seqs, species1)
                        with open(oma_seqs) as omaseq:
                            for line2 in omaseq:
                                if '>' in line2 and species1Id in line2:
                                    fnew.write(line2[:6] + '\n' + omaseq.next().replace('*', ''))
                                    break
    
                    if os.path.exists(oma_proteomes_dir + '/proteome_' + species2):
                        with open(oma_proteomes_dir + '/proteome_' + species2) as prot2:
                            for prot_line in prot2:
                                if '>' in prot_line and species2Id in prot_line:
                                    fnew.write(prot_line[:6] + '\n' + prot2.next().replace('*', ''))
                                    break
                    else:
                        print('Searching in file: ', oma_seqs, species2)
                        with open(oma_seqs) as omaseq:
                            for line2 in omaseq:
                                if '>' in line2 and species2Id in line2:
                                    fnew.write(line2[:6] + '\n' + omaseq.next().replace('*', ''))
                                    break
                # ProtTrace rather uses MAFFT linsi than muscle
                # To avoid adding another dependency, I uncommented the linsi command
                #os.system('muscle -quiet -in %s -out %s' %(filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)))
                os.system('{0} --quiet --phylipout --thread {1} {2} > {3}'.format(linsi, nr_processors, filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)))

                # Increments the species counter for the log
                c += 1
        # If the counter has never been incremented,
        # we can assume that the species pair is missing in the file
        if c == 1:
            print("ERROR: No orthologous pairs found between {0} and {1}!".format(species1,species2))

    # Collect all pairwise protein alignments, concatenate them
    # and calculate the summarized pairwise species distance
    print('Preprocessing complete..\nConcating the alignments..')
    concat_file = species1+'_'+species2+'.aln'
    # The concatAlignment script requires perl
    os.system('perl %s -in=%s -out=%s' %(concatAlignment, aln_dir, concat_file))
    print('Postprocessing concatenated alignment file..')
    post_process_concat_file = concat_file.replace('.aln', '.dup.aln')
    fnew = open(post_process_concat_file, 'w')
    with open(concat_file) as con:
        for line in con:
            if '>' in line:
                seq_id = line
                sequence = con.next()
                fnew.write(seq_id + sequence + seq_id.replace('>', '>dup_') + sequence)
    fnew.close()
    print('Converting alignment file to phylip format..')
    phyFile = post_process_concat_file.replace('.aln', '.phy')
#    os.system('%s -convert -output=phylip -infile=%s -outfile=%s' %(clustalw, post_process_concat_file, phyFile))
    
    print('Performing likelihood mapping..')
    with open('temp_puzzleParams.txt', 'w') as pp:
        pp.write(phyFile + '\nb\ne\nm\nm\nm\nm\nm\nm\ny\n')

    os.system('puzzle < %s' %puzzleParams)
    
    print('Parsing outfile and saving likelihood distance to')
    result_file = species1 + '_' + species2 + '.lik'
    with open(result_file,'w') as result:
        outdist = open(phyFile + '.dist').read().split('\n')
        result.write(outdist[3].split()[1])

    if target_dir is not None:
        # If the target_dir is just a directory, append the file name
        if add_filename == 1:
            target_dir_copy_path = os.path.join(os.path.abspath(target_dir), result_file)
        # If the target_dir looks like a filename on its own, use it directly
        else:
            target_dir_copy_path = target_dir
        with open(target_dir_copy_path,'w') as result:
            outdist = open(phyFile + '.dist').read().split('\n')
            result.write(outdist[3].split()[1])

    #if delete_temp: 
    # Deletes temporary files 
    #    if os.path.exists(result) and not len(open(result).read().split('\n')) == 0:
    #        os.system('rm -rf %s' %fa_dir)
    #        os.system('rm -rf %s' %aln_dir)
    #        os.system('rm %s*' %phyFile)
    #        os.system('rm -rf *.txt')
    #        os.system('rm -rf *.aln')

    os.chdir(root_dir)

    # END OF POSTPROCESS FUNCTION DEFINITION

# This defines the start of this script if someone wants to calculate
# distances separately with this script
if __name__ == "__main__":

    # Define and check for arguments
    try:
        opts, args = getopt.getopt(argv, "q:s:o:r:h", ["query=", "species=", "output=", "help"])
    except getopt.GetoptError:
        print('Invalid arguments:\nUsage:\tdistanceCalculation.py -q <query_species_id> -s <paired_species_id> -o <output_dir> [-help]')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h','--help'):
            print("USAGE:\tdistanceCalculation.py -q <query_species_id> -s <paired_species_id> -o <output_file>\n\t-q\t\tThe OMA ID of the query species\n\t-s\t\tThe OMA ID of the target species\n\t-o\t\tThe full output file path")
            sys.exit(2)
        elif opt in ('-q', '--query'):
            query = arg
        elif opt in ('-s','--species'):
            species = arg
        elif opt in ('-o','--output'):
            output = arg
        else:
            print('Invalid arguments:\nUsage:\tprotTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-help]')
            sys.exit(2)

    if len(sys.argv[1:]) == 0:
        print('ERROR:\tNo arguments entered for calculating species distances:\nUSAGE:\tdistanceCalculation.py -s <query species ID> -c <configFile> [-help]')
        sys.exit(2)

    # If this script is started on its own calculate the pairwise distance between the specified species pair
    calculate_protein_distances(query,species,None,output)
