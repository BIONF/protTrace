import os, sys
import glob
import random

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

    # DEBUG
    calculate_protein_distances(query,"NASVI",config,cache_dir,1,100)
    sys.exit()

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

def concatenate_alignment(species1, species2, alignment_count, alignment_directory, bootstrap_count=0):

    # The sequences of both species are separated for better handling
    # For bootstrapping, the first index tells us the aligned protein pair
    # The second index tells us the species
    sequences = ["",""]
    protein_pair = 0
    linecount = 0
    subalignment_positions = []
    subalignment_count = -1
    # For a sanity check, we compare the input species with the species
    # listed in the first two sequence lines of the file
    sanity_line_count = 0
    species_found = 0

    # Assuming PHYLIP format, collect all PHYLIP files in the directory
    print('Collecting all pairwise alignments for concatenation!')
    for current_alignment in range(1, alignment_count):
        with open("{0}/seq_{1}.aln".format(alignment_directory, str(current_alignment))) as phy_input:
            # The first line in the phylip file shows the first and
            # last original protein positions. We put these into a 
            # separate file for later phylogenomic diagnosis

            # The positions are converted from 1 based to 0 based counting
            # Hence the -1 substractions from the read-in numbers
            # Then, the first position of the next sequence needs to be
            # the last position of the previous one + 1
            if subalignment_count > -1:
                subalignment_begin = int(subalignment_positions[subalignment_count][1]) + 1
                subalignment_positions.append([subalignment_begin, subalignment_begin + int(phy_input.readline().strip().split()[1])])
            else:
               subalignment_positions.append([0, int(phy_input.readline().strip().split()[1]) - 1])
            subalignment_count += 1
            for line in phy_input:
                if line != "\n":
                    # The first two lines containing sequences
                    # are checked for whether their species IDs
                    # corresponds to the passed species IDs
                    if sanity_line_count < 2:
                        if line[0:11].strip() == species1:
                            species_found += 1
                        if line[0:11].strip() == species2:
                            species_found += 2
                        sanity_line_count += 1
                    if sanity_line_count == 3:
                        if species_found != 3:
                            print("ERROR: The passed species are not found in the alignment!\n{0}-{1}|{2}-{3}".format(species1,sequences[0][0:11].strip(),species2,sequences[1][0:11].strip()))
                        # Set the line count to anything that evades further operations
                        sanity_line_count = 4

                    # We add the species id into the indented area later
                    # Spaces must be stripped to generate a continuous 
                    # sequence that can be concatenated easily
                    sequences[linecount] += line[11:].strip().replace(" ","")
                    # The linecount is switched back and forth
                    if linecount == 0:
                        linecount = 1
                    else:
                        linecount = 0

    # Perform a bootstrap analysis on the alignment and note the created variance
    def create_bootstrap(protein_pairs, count):

        # Seed the singleton random module using the current system clock time
        # (by passing no parameter)
        seed = random.randrange(sys.maxsize)
        random.seed(seed)

        print("The seed for bootstrapping the alignment is: ", seed)
        with open('bootstrap_sample_rng_seed.R','w') as seed_file:
            seed_file.write('# The seed for sampling columns for pairwise distance\n# bootstrap analysis is a follows:\n# '+ str(seed) + '\n')

        # Generate a list of sequence indices to sample
        # The list is sampled with equal weights and replacement
        return [random.choices([i for i in range(protein_pairs)], k=protein_pairs) for c in range(count)]

    bootstrap_alignment_indices = create_bootstrap(len(sequences[0]),bootstrap_count)

    def sequence_pair_to_phylip(species1,species2,sequences,subalignment_positions,sampled_indices=None,bootstrap_count=None):

        # Here, we compile the bootstrapped sequences to continuous strings
        # The regular concatenated protein alignments are just made continuous
        # Both species remain separated
        local_sequences = ["",""]
        if sampled_indices is not None:
            if len(sampled_indices) != len(sequences[0]):
                print("ERROR: The bootstrap sequence has not the same length as the protein pair count!")
            for i in range(len(sampled_indices)):
                local_sequences[0] += sequences[0][sampled_indices[i]]
                local_sequences[1] += sequences[1][sampled_indices[i]]
        else:
            local_sequences = sequences

        def append_phylip_seq_with_spaces(begin,seq):
            # The sequence is indented to the right by 10 spaces
            appended_seq = " " * 10
            # The next 50 positions are inserted with one space between
            # 10 positions. The minimum function ensures that we do not
            # try to access empty sequence positions
            for i in range(begin,begin + min(len(seq) - begin,50),10):
                appended_seq += " " + seq[i:i+10]
                # Fill an open ending block with hyphens
                appended_seq += "-" * (10 - len(seq[i:i+10]))
            appended_seq += "\n"
            return appended_seq
    
        def append_phylip_seq_with_spaces_with_species(species,seq):
            appended_seq = species + " " * (10 - len(species))
            for i in range(0,50,10):
                appended_seq += " " + seq[i:i+10]
            appended_seq += "\n"
            return appended_seq
    
        # Build the PHYLIP file from the continuous alignment
        alignment = ""
        # The first line contains the entire length of the original protein sequence
        # Here, we just copy the first and last position
        alignment += "4 " + str(sum([sub[1] for sub in subalignment_positions])) + "\n"
        # The next 2 lines contain the respective species ids of each line
        # We duplicate the sequence to generate an imaginative 
        # 4 species alignment. We need this to calculate a 
        # tree and calculate the distance between sister 
        # clades (each clade is the same species twice)
        alignment += append_phylip_seq_with_spaces_with_species(species1,local_sequences[0])
        alignment += append_phylip_seq_with_spaces_with_species(species1 + "_dub",local_sequences[0])
        alignment += append_phylip_seq_with_spaces_with_species(species2,local_sequences[1])
        alignment += append_phylip_seq_with_spaces_with_species(species2 + "_dub",local_sequences[1])
        # Adds spacing between alignment blocks
        alignment += "\n"
        # This will fill the rest of the alignment file. We assume that both sequences
        # are equally long
        linecount = 0
        for i in range(50, len(local_sequences[0]), 50):
            alignment += append_phylip_seq_with_spaces(i,local_sequences[0])
            alignment += append_phylip_seq_with_spaces(i,local_sequences[0])
            alignment += append_phylip_seq_with_spaces(i,local_sequences[1])
            alignment += append_phylip_seq_with_spaces(i,local_sequences[1])
            # Adds spacing between alignment blocks
            alignment += "\n"
            if linecount == 0:
                linecount = 1
            else:
                linecount = 0

        # Assemble the name of the concatenation file
        # The bootstrap concatenation file is named as such
        concat_file_name = species1 + '_' + species2
        if sampled_indices is not None:
            concat_file_name += '_bootstrap_{0}'.format(bootstrap_count)
        concat_file_name += '.phy'
 
        with open(concat_file_name, 'w') as concat_file:
            concat_file.write(alignment)
        # Only the non-bootstrap result file receives a
        # subalignment_positions file
        if sampled_indices is None:
            with open('subalignment_positions.txt', 'w') as subpositions_file:
                for line in subalignment_positions:
                    # Columns are separated with spaces
                    subpositions_file.write(' '.join([str(l) for l in line]) + "\n")

    print('Concatenate and duplicate the pairwise aligned sequences!')
    # Generate the main concatenated and duplicated alignment for calculating
    # the pairwise species distance
    sequence_pair_to_phylip(species1,species2,sequences,subalignment_positions)

    print('Concatenate and duplicate the bootstrapped pairwise alignments!')
    # Generate the bootstrapped versions of the concatenated alignment for 
    # estimating the distance variance
    for i in range(len(bootstrap_alignment_indices)):
        sequence_pair_to_phylip(species1,species2,sequences,subalignment_positions,bootstrap_alignment_indices[i],i)
    
### Calculates the pairwise species maximum likelihood distance 
### between species1 and species2. The distance is copied to the
### target_dir. The config is loaded from ProtTrace's configure.py
### if this script is executed directly.
def calculate_protein_distances(species1,species2,config,target_dir,add_filename=1,bootstrap_count=0):

    print("Calculating the protein distance between {0} and {1}.".format(species1,species2))

    # Read the configuration of ProtTrace for paths
    oma_seqs = config.path_oma_seqs
    oma_pairs = config.path_oma_pairs
    concatAlignment = config.concat_alignments_script
    oma_proteomes_dir = config.path_distance_work_dir
    linsi = config.msa
    #clustalw = config.clustalw
    treepuzzle = "/share/applications/tree-puzzle/bin/puzzle"
    nr_processors = config.nr_processors
    delete_temp = config.delete_temp

    root_dir = os.getcwd()
    
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

    print(os.getcwd())
    
    # Align the set of pairwise oprthologs between the query
    # and the target species
    print('Preprocessing:\tParsing sequence pairs and aligning them...')
    prot_pairs = {}
    sequence_count = 1
    try:
        with open(oma_pairs,'r') as pair_mapping:
            for line in pair_mapping:
                # Identify the input species pair
                # among the orthologous protein pairs
                if species1 in line and species2 in line:
                    # Grab the protein ID
                    splitted_line = line.rstrip().split('\t')
                    if species1 in splitted_line[0]:
                        species_1_prot_id = splitted_line[0]
                        species_2_prot_id = splitted_line[1]
                    else:
                        species_2_prot_id = splitted_line[0]
                        species_1_prot_id = splitted_line[1]

                    # Store the protein IDs in a dictionary
                    # This way, we can immediately recognize the correct
                    # protein pair for one detected protein ID
                    # Each protein ID is a key to the same protein
                    # pair list
                    # The list also contains the sequence count to
                    # name each pairwise orthologous alignment file correctly
                    prot_pair = [species_1_prot_id, species_2_prot_id, sequence_count]
                    prot_pairs[species_1_prot_id] = prot_pair
                    prot_pairs[species_2_prot_id] = prot_pair

                    # Inform the user about the current progress
                    if sequence_count % 100 == 0:
                        print('Sequence Nr.:', sequence_count)
                                                      
                    #DEBUG
                    if sequence_count == 20:
                        break

                    sequence_count += 1

    except KeyboardInterrupt:
        sys.exit('The user interrupted the sequence search!')
    except FileNotFoundError:
        sys.exit('ERROR: The OMA pairs file is missing!') 

    # If the counter has never been incremented,
    # we can assume that the species pair is missing in the file
    if sequence_count == 1:
        print("ERROR: No orthologous pairs found between {0} and {1}!".format(species1,species2))

    ### Search the sequences of all protein pairs
    ### Each pair will be aligned directly after

    # All proteomes can be located in separate directories
    # that follow a systematic nomenclature
    # Otherwise, we assume that all proteins of the 
    # species can be found within the oma_seqs.fa file
    sequence_sources = set()
    if os.path.exists(oma_proteomes_dir + '/proteome_' + species1):
        sequence_sources.add(oma_proteomes_dir + '/proteome_' + species1)
    else:
        sequence_sources.add(oma_seqs)
    if os.path.exists(oma_proteomes_dir + '/proteome_' + species2):
        sequence_sources.add(oma_proteomes_dir + '/proteome_' + species2)
    else:
        sequence_sources.add(oma_seqs)

    # We iterate through both sequence sources. Since we
    # iterate a set, we save iterations if both point to
    # the same file
    try:
        for sequence_source in sequence_sources:
            with open(sequence_source, 'r') as sequence_source_file:
                for line in sequence_source_file:
                    if '>' in line:
                        line_id = line.replace('>', '')
                        if line_id in prot_pairs:
                            # The third element in each prot_pair list is their numeric sequence ID
                            # In append mode, the file is created if it does not exist yet
                            with open(fa_dir + '/seq_' + str(prot_pairs[line_id][2]) + '.fa', 'a') as alignment_precurser:
                                alignment_precurser.write(line[:6] + '\n' + sequence_source_file.readline().replace('*', ''))
    except FileNotFoundError:
        print('ERROR: The FASTA file that contains the sequences is missing!')

    # We perform a global alignment of each orthologous pair
    # For this step, we only need to know the number of pairs
    # to access the fasta files in the fasta directory
    # ProtTrace rather uses MAFFT linsi than muscle
    # To avoid adding another dependency, I uncommented the linsi command
    #os.system('muscle -quiet -in %s -out %s' %(filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)))
    for pair in range(1, sequence_count):
        filename = fa_dir + '/seq_' + str(pair) + '.fa'
        os.system('{0} --quiet --phylipout --thread {1} {2} > {3}'.format(linsi, nr_processors, filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)))

    # Collect all pairwise protein alignments, concatenate them
    # and calculate the summarized pairwise species distance
    print('Preprocessing complete..\nConcating the alignments..')

    # Concatenate and bootstrap the pairwise protein alignments
    concatenate_alignment(species1, species2, sequence_count, aln_dir, bootstrap_count)

    # Anticipate the output file names of the concatenated and
    # the bootstrapped concatenated alignments
    concat_files = [species1+'_'+species2+'.phy']
    for i in range(bootstrap_count):
        concat_files.append(species1+'_'+species2+'_bootstrap_{0}.phy'.format(i))

    # The concatAlignment script requires perl
#    os.system('perl %s -in=%s -out=%s' %(concatAlignment, aln_dir, concat_file))
    print('Postprocessing the concatenated alignment file..')

    # Executes TreePUZZLE to calculate the tree distance between the species pair
    def calculate_pairwise_distance(concat_file,treepuzzle):
        
        print('Performing likelihood mapping..')
        # Prepare the parameter file for TreePUZZLE
        with open('temp_puzzleParams.txt', 'w') as pp:
            #pp.write(concat_file + '\nb\ne\nm\nm\nm\nm\nm\nm\ny\n')
            pp.write(concat_file + '\ne\nm\nm\nm\nm\nm\nm\ny\n')
   
        # Execute TreePUZZLE
        os.system('{0} < temp_puzzleParams.txt >/dev/null'.format(treepuzzle))

    # Record the calculated main and bootstrap distances
    distances = []

    for concat_filename in concat_files:
        calculate_pairwise_distance(concat_filename,treepuzzle)
        # Read the calculated distance into a list of distances
        with open(concat_filename + '.dist','r') as concat:
            next(concat)
            distances.append(concat.readline().split()[3])

    # Write the main computed pairwise species distance to a result file
    result_file = concat_files[0].replace(".phy",".lik")
    with open(result_file,'w') as result:
        result.write(distances[0] + '\n')

    # Copy the main output file to the given target directory, if given
    # The target directory is originally designed to be the ProtTrace cache directory
    if target_dir is not None:
        # If the target_dir is just a directory, append the file name
        if add_filename == 1:
            target_dir_copy_path = os.path.join(os.path.abspath(target_dir), concat_files[0].replace(".phy",".lik"))
        # If the target_dir looks like a filename on its own, use it directly
        else:
            target_dir_copy_path = target_dir
        with open(target_dir_copy_path,'w') as result:
            result.write(distances[0] + '\n')

    # Write the main distance and its bootstrap distances into a table
    # for followup statistical analyses
    # The bootstrap column is a R friendly boolean to separate the main distance
    # from the bootstrap values
    separator = '\t'
    columns = [separator.join(['Species_1','Species_2','Distance','Bootstrap'])]
    for d in range(len(distances)):
        columns.append(separator.join([species1,species2,distances[d],'FALSE' if d == 0 else 'TRUE']))
    with open(concat_files[0].replace('.phy', '_bootstrap.lik'),'w') as analysis_output:
        analysis_output.write('\n'.join(columns) + '\n')

    if delete_temp: 
        # Deletes temporary files 
        if os.path.exists(result_file) and not len(open(result_file).read().split('\n')) == 0:
            os.system('rm -rf {0} {1} temp_puzzleParams.txt'.format(fa_dir, aln_dir))
            os.system('rm -rf {0}_{1}.phy*'.format(species1, species2))
            # If there are more than one alignment, we assume these to be bootstraps
            # This check rather exists to prevent file not found errors
            if len(concat_files) > 1:
                os.system('rm -rf {0}_{1}_bootstrap_*'.format(species1, species2))


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