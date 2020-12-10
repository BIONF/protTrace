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
    calculate_protein_distances(query,"NASVI",config,cache_dir,1,10)
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

    # Assuming PHYLIP format, append all PHYLIP files in the directory
    for current_alignment in range(1, alignment_count):
        with open("{0}/seq_{1}.aln".format(alignment_directory, str(current_alignment))) as phy_input:
            # The first line in the phylip file shows the first and
            # last original protein positions. We put these into a 
            # separate file for later phylogenomic diagnosis
            if subalignment_count > -1:
                subalignment_positions.append(subalignment_positions[subalignment_count].split()[1] + " " + phy_input.readline().strip().split()[1])
            else:
               subalignment_positions.append("0" + " " + phy_input.readline().strip().split()[1])
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
        alignment += "4 " + subalignment_positions[-1].split()[1] + "\n"
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
        concat_file_name = species1 + '_' + species2
        if sampled_indices is not None:
            concat_file_name += '_bootstrap_{0}'.format(bootstrap_count)
        concat_file_name += '.phy'
 
        with open(concat_file_name, 'w') as concat_file:
            concat_file.write(alignment)
        if sampled_indices is None:
            with open('subalignment_positions.txt', 'w') as subpositions_file:
                for line in subalignment_positions:
                    subpositions_file.write(line + "\n")

    # Generate the main concatenated and duplicated alignment for calculating
    # the pairwise species distance
    sequence_pair_to_phylip(species1,species2,sequences,subalignment_positions)
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
    sequence_count = 1
    with open(oma_pairs,'r') as op:
        for line in op:
            if species1 in line and species2 in line:
                print('Sequence Nr.:', sequence_count)

                # Grab the protein ID
                splitted_line = line.rstrip().split('\t')
                if species1 in splitted_line[0]:
                    species1Id = splitted_line[0]
                    species2Id = splitted_line[1]
                else:
                    species2Id = splitted_line[0]
                    species1Id = splitted_line[1]

                # Search for both protein sequences
                filename = fa_dir + '/seq_'+str(sequence_count)+'.fa'
                with open(filename, 'w') as alignment_precurser:
                    if os.path.exists(oma_proteomes_dir + '/proteome_' + species1):
                        with open(oma_proteomes_dir + '/proteome_' + species1,'r') as prot1:
                            for prot_line in prot1:
                                if '>' in prot_line and species1Id in prot_line:
                                    alignment_precurser.write(prot_line[:6] + '\n' + prot1.readline().replace('*', ''))
                                    break
    
                    else:
                        print('Searching in file: ', oma_seqs, species1)
                        with open(oma_seqs,'r') as omaseq:
                            for line2 in omaseq:
                                if '>' in line2 and species1Id in line2:
                                    alignment_precurser.write(line2[:6] + '\n' + omaseq.readline().replace('*', ''))
                                    break
    
                    if os.path.exists(oma_proteomes_dir + '/proteome_' + species2):
                        with open(oma_proteomes_dir + '/proteome_' + species2,'r') as prot2:
                            for prot_line in prot2:
                                if '>' in prot_line and species2Id in prot_line:
                                    alignment_precurser.write(prot_line[:6] + '\n' + prot2.readline().replace('*', ''))
                                    break
                    else:
                        print('Searching in file: ', oma_seqs, species2)
                        with open(oma_seqs,'r') as omaseq:
                            for line2 in omaseq:
                                if '>' in line2 and species2Id in line2:
                                    alignment_precurser.write(line2[:6] + '\n' + omaseq.readline().replace('*', ''))
                                    break
                # ProtTrace rather uses MAFFT linsi than muscle
                # To avoid adding another dependency, I uncommented the linsi command
                #os.system('muscle -quiet -in %s -out %s' %(filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)))
                os.system('{0} --quiet --phylipout --thread {1} {2} > {3}'.format(linsi, nr_processors, filename, filename.replace('.fa', '.aln').replace(fa_dir, aln_dir)))

                # Increments the sequence counter
                sequence_count += 1

                #DEBUG
                if sequence_count == 8:
                    break

        # If the counter has never been incremented,
        # we can assume that the species pair is missing in the file
        if sequence_count == 1:
            print("ERROR: No orthologous pairs found between {0} and {1}!".format(species1,species2))

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
    print('Postprocessing concatenated alignment file..')
    
    # Duplicate both sequences to 4 to make a minimum tree possible
#    post_process_concat_file = concat_file.replace('.aln', '.dup.aln')
#    fnew = open(post_process_concat_file, 'w')
#    with open(concat_file) as con:
#        for line in con:
#            if '>' in line:
#                seq_id = line
#                sequence = con.readline()
#                fnew.write(seq_id + sequence + seq_id.replace('>', '>dup_') + sequence)
#    fnew.close()

    print('Converting alignment file to phylip format..')
    #phyFile = post_process_concat_file.replace('.aln', '.phy')
#    os.system('%s -convert -output=phylip -infile=%s -outfile=%s' %(clustalw, post_process_concat_file, phyFile))

    # Executes TreePUZZLE to calculate the tree distance between the species pair
    def calculate_pairwise_distance(concat_file,treepuzzle):
        
        print('Performing likelihood mapping..')
        # Prepare the parameter file for TreePUZZLE
        with open('temp_puzzleParams.txt', 'w') as pp:
            pp.write(concat_file + '\nb\ne\nm\nm\nm\nm\nm\nm\ny\n')
    
        # Execute TreePUZZLE
        os.system('{0} < temp_puzzleParams.txt >/dev/null'.format(treepuzzle))

    for concat_filename in concat_files:

        calculate_pairwise_distance(concat_filename,treepuzzle)

        # Write the distance from TreePuizzle's output into a simpler format
        result_file = concat_filename.replace(".phy",".lik")
        with open(result_file,'w') as result:
            with open(concat_filename + '.dist') as concat:
                # We rather take line 1 (counted from 0), column 4, because
                # we reach it faster than line 4, column 1
                next(concat)
                result.write(concat.readline().split()[4] + "\n")

    # Only the main concatenated alignment is taken as the average distance
    # Bootstraps are analyzed somewhere else
#    print('Parsing outfile and saving likelihood distance to')
#    result_file = species1 + '_' + species2 + '.lik'
    
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
            with open(concat_files[0] + '.dist') as concat:                          
                # We rather take line 1 (counted from 0), column 4, because
                # we reach it faster than line 4, column 1
                next(concat)
                result.write(concat.readline().split()[4] + "\n")

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
