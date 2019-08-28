import os, sys
import glob
import subprocess
from multiprocessing import Pool

# Script to run HAMSTR search to add orthologs to the core-orthologs groups derived from OMA database
# Input: 1.Path to HAMSTR folder
# 	 2.Path to orthologs group fasta file
#	 3.Protein id
#	 4.Mapping file -> HAMSTR id to OMA id
# Output: Adds orthologs found to the core-ortholog file
# NOTE: The $path and $output variables in the hamstr.pl file has been edited.

# function to remove if pre-existing run data for a protein id is present
def remove_old_dir():
    if os.path.exists(hamstr + '/core_orthologs/' + protein_id):
        os.system('rm -Rf %s/core_orthologs/%s' %(hamstr, protein_id))
    if os.path.exists(hamstr + '/output/' + protein_id):
        os.system('rm -Rf %s/output/%s' %(hamstr, protein_id))

# function to make new directories for HAMSTR search
def make_new_dir():

    os.system('mkdir %s' %core_ortholog_prot_dir)
    os.system('mkdir %s/aln_dir' %core_ortholog_prot_dir)
    os.system('mkdir %s/hmm_dir' %core_ortholog_prot_dir)
    os.system('mkdir %s' %output)

# Function to read the core orthologs file and make a copy of it in the HAMSTR folder
# Edits the fasta file how in format required for the HAMSTR run
def copy_edit_core_ortholog():
    global fasta_file
    fasta_file = core_ortholog_prot_dir + '/' + protein_id + '.fa'
    s1 = open(core_ortholog).read().split('\n')
    fnew = open(fasta_file, 'w')
    flag = True ###	Check whether HaMStR has a starting sequence	###
    for i in range(0, len(s1)-1, 2):
        if s1[i][0] == '>':
            sequence = s1[i+1].replace('*', '').replace('X', '')
            s2 = open(hamstr_map_oma).read().split('\n')
            for j in range(len(s2)-1):
                if s2[j].split()[-1] == s1[i][1:6]:
                    hamstrProtId = s2[j].split()[0]
                    print (hamstrProtId + '\n')
                    global blastDir
                    if not hamstr_env == "":
                        blastDir = hamstr + '/blast_dir_' + hamstr_env
                    else:
                        blastDir = hamstr + '/blast_dir'
                    for dirs in glob.glob(blastDir + '/*'):
                        if hamstrProtId in dirs:
                            blastFastaFile = dirs + '/' + dirs.split('/')[-1] + '.fa'
                            break
                    s3 = open(blastFastaFile).read().split('\n')
                    runBlast = True # If the sequence is not found in our database then, take the top blast hit #
                    for k in range(len(s3) - 1):
                        if sequence == s3[k].replace('*', '').replace('X', ''):
                            fnew.write(s1[i][0] + protein_id + '|' + dirs.split('/')[-1] + '|' + s3[k-1][1:] + '\n' + sequence + '\n')
                            flag = False
                            runBlast = False
                            break

                    if runBlast:
                        print 'No matching sequence found in hamstr blast directory! Running BLAST search now..'
                        currentWorkDir = os.getcwd()
                        print (currentWorkDir + '\n')
                        tempDir = 'temp_blast_' + hamstrProtId
                        if not os.path.exists(tempDir):
                            os.mkdir(tempDir)
                        os.chdir(tempDir)
                        temp_query = open('temp_query.fa', 'w')
                        temp_query.write(s1[i] + '\n' + s1[i+1])
                        temp_query.close()
                        print (blastFastaFile + 'blastFastafile \n')
                        os.system('ln -s %s temp_proteome.fa' %(blastFastaFile))
                        os.system('%s -i temp_proteome.fa' %(format_db))
                        os.system('%s -query temp_query.fa -db temp_proteome.fa -evalue 0.00001 -outfmt 6 -max_target_seqs 1 -out temp_out.txt' %(blastp))

                        if os.path.exists('temp_out.txt') and len(open('temp_out.txt').read().split('\n')) > 1:
                            hit_id = open('temp_out.txt').read().split('\n')[0].split('\t')[1]
                            fnew.write(s1[i][0] + protein_id + '|' + dirs.split('/')[-1] + '|' + hit_id + '\n' + sequence + '\n')
                        os.chdir(currentWorkDir)

                        if delTemp:
                            os.system('rm -rf %s' %tempDir)
                    break
    fnew.close()

    if flag:
        return False
    else:
        return True

# Function to run mafft and hmmbuild on the edited core ortholog fasta file
def align_hmm_run():
    global aln_file, hmm_file
    aln_file = core_ortholog_prot_dir + '/aln_dir/' + protein_id + '.aln'
    hmm_file = core_ortholog_prot_dir + '/hmm_dir/' + protein_id + '.hmm'

    f = open(fasta_file).read().split('\n')
    if len(f) > 3:
        os.system('linsi %s > %s' %(fasta_file, aln_file))
        os.system('hmmbuild %s %s' %(hmm_file, aln_file))
    else:
        os.system('hmmbuild %s %s' %(hmm_file, fasta_file))

def actual_hamstr_computation(dirs):
    if not ".fa" in dirs and not ".sql" in dirs:
        hamstr_name = dirs.split('/')[-1]
        if hamstr_name in hamstr_species_in_tree:
            file_name = dirs + '/' + dirs.split('/')[-1] + '.fa'
            if os.path.exists(file_name + '.mod'):
                file_name = file_name + '.mod'

            if not include_paralogs:
                command = '%s/bin/hamstr.pl -central -sequence_file=%s -taxon=misc -hmmset=%s -strict -checkCoorthologsRef -representative -outpath=%s -scoreThreshold -blastpath=%s' %(hamstr, file_name, protein_id, output, blastDir)
            if include_paralogs:
                command = '%s/bin/hamstr.pl -central -sequence_file=%s -taxon=misc -hmmset=%s -strict -checkCoorthologsRef -outpath=%s -coreThreshold -blastpath=%s' %(hamstr, file_name, protein_id, output, blastDir)
            #command = '%s/bin/hamstr.pl -central -sequence_file=%s -taxon=misc -hmmset=%s -strict -representative -outpath=%s -scoreThreshold -blastpath=%s' %(hamstr, file_name, protein_id, output, blastDir)
            print 'HaMStR run with command: ', command
            try:
                subprocess.call(command, shell=True)
            except KeyboardInterrupt as e:
                raise Exception
            except:
                raise Exception

# Function to run HAMSTR
def hamstr_run():
    global hamstr_species_in_tree
    hamstr_species_in_tree = []
    s2 = open(hamstr_map_oma).read().split('\n')
    for i in range(len(s2) - 1):
        hamstr_species_in_tree.append(s2[i].split('\t')[0])

    if not hamstr_env == "":
        genomeDirRun = '/genome_dir_%s/*' %hamstr_env
    else:
        genomeDirRun = '/genome_dir/*'

    dirsList = []
    for dirs in glob.glob(hamstr + genomeDirRun):
        dirsList.append(dirs)

    try:
        pool = Pool(processes=nr_proc)
        pool.map(actual_hamstr_computation, dirsList)
    except KeyboardInterrupt as e:
        pool.terminate()
        pool.join()
        sys.exit(e)
    except:
        print("ERROR: Multiprocessing step <-> HaMStR search.")
        pass


# Function to add orthologs to the core orthologs set
def write_output():
    sequenceDict = {} #This dictionary helps to avoid duplicated sequences in final orthologs file.

    #print("HaMStR Output dir: ", output)
    count = 0
    for f2 in glob.glob(output+'/*.out'):
        count += 1

    #print("HaMStR .out count: ", count)

    if count > 0:
        speciesList = {} #Variable to store the name of the core-orthologs species (eg -> YEAST00011 would be stored in the list as 'YEAST')

        s1 = open(core_ortholog).read().split('\n')
        fnew = open(core_ortholog, 'w')

        s4 = open(hamstr_map_oma).read().split('\n')
        present_species = []

        for i in range(len(s4) - 1):
            present_species.append(s4[i].split()[0])

        # Re-writing the original core orthologs set
        for i in range(len(s1)-1):
            if '>' in s1[i]:
                target_species_id = s1[i].split('>')[-1]
                targetSequence = s1[i+1].replace('X', '').replace('*', '')
                if not target_species_id in speciesList.keys():
                    speciesList[target_species_id] = []
                    speciesList[target_species_id].append(targetSequence)
                    fnew.write(s1[i].replace('>', '>1_') + '\n' + targetSequence + '\n')
                else:
                    if not targetSequence in speciesList[target_species_id]:
                        speciesList[target_species_id].append(targetSequence)
                        fnew.write(s1[i].replace('>', '>1_') + '\n' + targetSequence + '\n')
        #print co
        # Reading the .out files and adding the sequences to the core orthologs set
        for files in glob.glob(output+'/*.out'):
            #print 'Writing file: ', files
            s2 = open(files).read().split('\n')

            if '_' in files:
                #hamstr_name = '_'.join(files.split('/')[-1].split('.strict.out')[0].split('_')[1:])
                hamstr_name = files.split('/')[-1].split('_')[1]
            else:
                #hamstr_name = '_'.join(files.split('/')[-1].split('.out')[0].split('_')[1:])
                print 'ERROR: Please check .out files produced by HaMStR run ("_" MISSING in .out file names)'
                hamstr_name = files.split('/')[-1]

            s3 = open(hamstr_map_oma).read().split('\n')

            oma_name = ''
            for j in range(len(s3)-1):
                if s3[j].split()[0] == hamstr_name:
                    #oma_name = s3[j].split()[3] + s2[0].split('|')[3][:5]
                    oma_name = s3[j].split()[-1]
                    break
            #print oma_name[:5]
            if not oma_name in speciesList.keys() and hamstr_name in present_species:
                speciesList[oma_name] = []
                c = 1
                for i in range(len(s2) - 1):
                    targetSequence = s2[i].split('|')[-1].replace('X', '').replace('*', '')
                    if not targetSequence in speciesList[oma_name]:
                        speciesList[oma_name].append(targetSequence)
                        fnew.write('>' + str(c) + '_' + oma_name + '\n' + targetSequence + '\n')
                        c += 1
            elif oma_name in speciesList.keys() and hamstr_name in present_species:
                c = len(speciesList[oma_name]) + 1
                for i in range(len(s2) - 1):
                    targetSequence = s2[i].split('|')[-1].replace('X', '').replace('*', '')
                    #print speciesList[oma_name]
                    #print s2[i]
                    if not targetSequence in speciesList[oma_name]:
                        speciesList[oma_name].append(targetSequence)
                        fnew.write('>' + str(c) + '_' + oma_name + '\n' + targetSequence + '\n')
                        c += 1
                    #print speciesList[oma_name]
        fnew.close()
    else:
        print 'No orthologs to be added!!'

# Edit the orthologs file to just contain species which are present in the species tree
def finalEditOrthologsSeqs():
    s1 = open(core_ortholog).read().split('\n')
    fnew = open(core_ortholog, 'w')

    species_in_tree = []
    s2 = open(hamstr_map_oma).read().split('\n')
    for i in range(len(s2) - 1):
        species_in_tree.append(s2[i].split('\t')[-1])

    for i in range(0, len(s1) - 2, 2):
        if s1[i][0] == '>':
            if '_' in s1[i]:
                if s1[i][3:8] in species_in_tree:
                    fnew.write(s1[i] + '\n' + s1[i+1] + '\n')
            else:
                if s1[i][1:6] in species_in_tree:
                    fnew.write(s1[i] + '\n' + s1[i+1] + '\n')
    fnew.close()

# Remove the HaMStR created directories
def removeHaMStRdirs():
    if os.path.exists(hamstr + '/core_orthologs/' + protein_id):
        os.system('rm -Rf %s/core_orthologs/%s' %(hamstr, protein_id))
    if os.path.exists(hamstr + '/output/' + protein_id):
        os.system('rm -Rf %s/output/%s' %(hamstr, protein_id))

# The main method to run HaMStR search
def main(hamstrFile, coreOrtholog, protId, hamstrMapOma, makeblastdb, blast, del_temp, hamstr_environment, include_para, nr_processors):
    print '##### Running HaMStR search #####'

    global hamstr, core_ortholog, protein_id, hamstr_map_oma, core_ortholog_prot_dir, output, fasta_file, aln_file, hmm_file, format_db, blastp
    global core_ortholog_prot_dir, output, delTemp, hamstr_env, include_paralogs, nr_proc
    protein_id = protId
    core_ortholog = coreOrtholog
    hamstr = hamstrFile
    core_ortholog_prot_dir = hamstr + '/core_orthologs/' + protein_id
    output = hamstr + '/output/' + protein_id
    hamstr_map_oma = hamstrMapOma
    format_db = makeblastdb
    blastp = blast
    delTemp = del_temp
    hamstr_env = hamstr_environment
    include_paralogs = include_para
    nr_proc = nr_processors

    print 'HaMStR processing Step 1: Removing old directories'
    remove_old_dir()
    print 'HaMStR processing Step 2: Making new directories'
    make_new_dir()
    print 'HaMStR processing Step 3: Edit core orthologs file'
    success = copy_edit_core_ortholog()
    #success = True
    print 'HaMStR processing Step 4: Creating alignment and HMMER file'
    if success:
        align_hmm_run()
        print 'HaMStR processing Step 5: Running HAMSTR'
        try:
            #pass
            hamstr_run()
        except KeyboardInterrupt:
            sys.exit('Keyboard interruption by user!!!')
        print 'HaMStR processing Step 6: Adding orthologs to the core set'
        write_output()
        print 'HaMStR processing Step 7: Editing the orthologs and keeping only those present in the species tree'
        finalEditOrthologsSeqs()
    print 'HaMStR processing Step 8: Removing the directories created by HaMStR'
    if delTemp:
        removeHaMStRdirs()

    return success






