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
# (cut down every fasta header to 5 character long names)
#def rename_orth_file():
#   fnew = open('temp_orth_%s.fa' %protein_id, 'w')
#   for i in range(0, len(orth_file) - 1, 2):
#       fnew.write(orth_file[i][:6] + '\n' + orth_file[i+1] + '\n')
#   fnew.close()

# Perform MSA of the new renamed ortholog sequences file (MAFFT linsi)
def msa_convert():
    if os.path.exists(phy_file):
        print("Reusing existing alignment file: phy_file")
    else:
        print("Generating MSA")
        #subprocess.run([linsi,"--phylipout","ogSeqs_{0}.fa".format(protein_id),"ogSeqs_{0}.phy".format(protein_id)])
        os.system('{0} --phylipout ogSeqs_{1}.fa > ogSeqs_{3}.phy'.format(linsi, protein_id, protein_id))

########### added by ingo to get rid of raxml dependency
def run_iqtree():
    if makeTree:
        if reuse_cache and os.path.exists('ogSeqs_' + protein_id + '.phy.treefile') and os.path.exists('ogSeqs_' + protein_id + '.phy.ckp.gz'):
            print('ML tree already exists. Reusing it.')
        else:
            os.system('rm -rf RAxML_*')
            os.system('iqtree -nt %s -s ogSeqs_%s.phy -m %s -keep-ident -redo' %(nr_proc, protein_id, aaMatrix))
        #print('complete..')
# Remove all the temp files generated
def rm_temp():

    os.system('rm temp_parameters_%s.txt' %protein_id)
    os.system('rm maxLikDist_%s.txt' %protein_id)

# Calculate the scaling factor based on maximum likelihood distances
def scalingFactorMax():
    scales = []
    # Generate maximum likelihood distance file for orthologs
    ### NOTE THE FILE USED HERE!!!!*******************************************!!!!
    outfile = open('ogSeqs_%s.phy.mldist' %protein_id).read().split('\n')
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
                species1 = line

            ###
            ### This block is needed when max likelihood matrix do not have OMA identifiers ###
            ###
            #print(species1)
            #for j in range(len(hamstrFile) - 1):
            #	if species1 == hamstrFile[j].split('\t')[3]:
            #		hamstr1 = hamstrFile[j].split('\t')[0]
            #		break

            for k in range(i + 1, len(orthMaxFile) - 1):
                line = orthMaxFile[k].split('\t')[0]
                if '_' in line:
                    species2 = line.split('_')[1]
                else:
                    species2 = line

                #for j in range(len(hamstrFile) - 1):
                #	if species2 == hamstrFile[j].split('\t')[3]:
                #		hamstr2 = hamstrFile[j].split('\t')[0]
                #		break
                #print(species1, species2)
                maxDistOrth = float(orthMaxFile[i].split('\t')[k + 1])
                maxDistSpecies = 0 #Maximum likelihood distance from species max. likelihood matrix
                #print(maxDistOrth)
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
                    #print('yes')
                    # Checking for the likelihood score in cache directory
                    if os.path.exists(cacheDir + '/' + species1 + '_' + species2 + '.lik'):
                        maxDistSpecies = float(open(cacheDir + '/' + species1 + '_' + species2 + '.lik').read().split('\n')[0])
                    elif os.path.exists(cacheDir + '/' + species2 + '_' + species1 + '.lik'):
                        maxDistSpecies = float(open(cacheDir + '/' + species2 + '_' + species1 + '.lik').read().split('\n')[0])
                    else:
                        #print('No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2))
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
                            #print('No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2))
                            mlPresent = False
                            #maxDistSpecies = 1.00

                #print(maxDistSpecies)
                if mlPresent and not maxDistSpecies == 0:
                    scales.append(maxDistOrth / maxDistSpecies)
                elif mlPresent and maxDistSpecies == 0:
                    #scales.append(1.00)
                    pass
    except:
        print('### ERROR: Scaling factor calculation had an error ###')
        sys.exit('Maximum likelihood files are invalid!')

    if len(scales) >= 1:
        return median(scales)
    else:
        return sf

# Main module for running tree reconstruction
def main(Linsi, Orthologs, AaMatrix, Protein_id, Map_file, Species_maxLikMatrix, Scale_file, Tree_file, delTemp, defScale, cache_dir, ortholog_tree_reconstruction, nr_processors, use_cache):

    global linsi, orth_file, aaMatrix, protein_id, map_file, species_maxLikMatrix, scaleFile, treeFile, phy_file, sf, cacheDir, makeTree, nr_proc, reuse_cache
    cacheDir = cache_dir
    linsi = Linsi
    aaMatrix = AaMatrix
    protein_id = Protein_id
    map_file = Map_file
    species_maxLikMatrix = Species_maxLikMatrix
    orth_file = open(Orthologs).read().split('\n')
    scaleFile = Scale_file
    treeFile = Tree_file
    phy_file = 'ogSeqs_' + protein_id + '.phy'
    sf = float(defScale)
    makeTree = ortholog_tree_reconstruction
    nr_proc = nr_processors
    reuse_cache = use_cache

    #sf = 1.00
    if len(orth_file) > 7:
        try:
#           rename_orth_file()
            msa_convert()
            run_iqtree()
            if reuse_cache and os.path.exists(scaleFile):
                print('Pre-computed scaling factor found. Reusing it.')
            else:
                sf = scalingFactorMax()
                print('Scaling factor: {0}'.format(sf))
            if delTemp:
                rm_temp()
        except:
            print('### ERROR: Some step in the tree reconstruction was invalid!! ###')
            pass
    else:
        print('Using default scaling factor: {0}'.format(sf))

    with open(scaleFile, 'w') as fnew:
        fnew.write(str(sf))
