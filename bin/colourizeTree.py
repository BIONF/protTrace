import os, sys
import math
import subprocess

### Script to colourize the species according to the traceability values
### Reads in the generated nexus file by PTP and edits it according to our need
### Represents the traceabilities as different colours on the leaves / tip labels
### INPUT: 1. Nexus tree file
###	   2. HaMStR mapping file
###	   3. Species Id
###	   4. Species tree
###	   5. Decay result file

def calculateMaxLikDist(species1, species2):
    #print(species1, species2)
    speciesMaxFile = open(sp_max_file).read().split('\n')
    for j in range(len(hamstrMapFile) - 1):
        if species1 == hamstrMapFile[j].split('\t')[1]:
            #hamstr1 = hamstrMapFile[j].split('\t')[0]
            oma1 = hamstrMapFile[j].split('\t')[-1]
            break
    for j in range(len(hamstrMapFile) - 1):
        if species2 == hamstrMapFile[j].split('\t')[1]:
            #hamstr2 = hamstrMapFile[j].split('\t')[0]
            #kingdom = hamstrMapFile[j].split('\t')[-1]
            oma2 = hamstrMapFile[j].split('\t')[-1]
            break
    #print(hamstr1, hamstr2)
    #print(oma1, oma2)
    flag1 = True #To check if the species are present in the species ML dist matrix file. Otherwise, parse likelihood from cache directory.
    flag2 = True

    #*** CHECK HERE FOR THE SEPARATION IDENTIFIER IN SPECIES LIKELIHOOD MATRIX FILE ('\t' or '|') ***#
    sep = '\t'

    for l in range(len(speciesMaxFile) - 1):
        if speciesMaxFile[l].split(sep)[0] == oma1:
            rowIndex = l
            flag1 = False
        elif speciesMaxFile[l].split(sep)[0] == oma2:
            columnIndex = l
            flag2 = False

    if flag1 or flag2:
        print('Checkpoint 2 crossed')
        # Checking for the likelihood score in cache directory
        if os.path.exists(cacheDir + '/' + oma1 + '_' + oma2 + '.lik'):
            return float(open(cacheDir + '/' + oma1 + '_' + oma2 + '.lik').read().split('\n')[0])
        elif os.path.exists(cacheDir + '/' + oma2 + '_' + oma1 + '.lik'):
            return float(open(cacheDir + '/' + oma2 + '_' + oma1 + '.lik').read().split('\n')[0])
        else:
            print('No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2))
            return 1.00
    else:
        if not speciesMaxFile[rowIndex].split(sep)[columnIndex] == "NA":
            return float(speciesMaxFile[rowIndex].split(sep)[columnIndex])
        else:
            print('No likelihood distance found between species: %s and %s. Using default likelihood distance of 1.0!' %(species1, species2))
            return 1.00

def colourize(speciesName, speciesId, nexusTreeFile, speciesTree, decayRate, decayPop, traceResults, matrixDict):
    tree = open(nexusTreeFile).read().split('\n')
    for i in range(len(tree) - 1):
        if tree[i] == "\ttaxlabels":
            startTaxa = i + 1
        if tree[i] == ";":
            stopTaxa = i
        if tree[i] == "begin trees;":
            treeLine = i + 2
            break
    fnew = open(nexusTreeFile.replace('.nexus', '_edit.nexus'), 'w')
    for i in range(startTaxa):
        fnew.write(tree[i] + '\n')

    for i in range(startTaxa, stopTaxa):
        tempSpecies = tree[i].split()[-1]
        if not tempSpecies == speciesName:
            try:
                colourCode, traceValue = getColourCode(speciesName, tempSpecies, decayRate, decayPop)
                omaName = "NA"

                for j in range(len(hamstrMapFile) - 1):
                    if tempSpecies == hamstrMapFile[j].split('\t')[1]:
                        omaName = hamstrMapFile[j].split('\t')[-1]
                        break
                for elements in matrixDict[omaName]:
                    newElement = elements + '#' + str(traceValue)
                    matrixDict[omaName].remove(elements)
                    matrixDict[omaName].insert(0, newElement)
                traceResults.write(speciesName + '\t' + tempSpecies + '\t' + str(traceValue) +'\n')
                fnew.write(tree[i] + '[&!color=#-' + colourCode + ']' + '\n')
            except:
                sys.exit("ERROR: Check species %s in species tree and used mapping files" %tempSpecies)
                pass

        else:
            fnew.write(tree[i] + '\n')
            for elements in matrixDict[speciesId]:
                newElement = elements + "#1"
                matrixDict[speciesId].remove(elements)
                matrixDict[speciesId].insert(0,newElement)

    for i in range(stopTaxa, treeLine):
        fnew.write(tree[i] + '\n')
    fnew.write(speciesTree.replace('\n', '') + '\n')
    for i in range(treeLine + 1, len(tree)):
        if tree[i] == '	set branchLabels.isShown=true;':
            fnew.write("	set branchLabels.isShown=false;" + '\n')
        else:
            fnew.write(tree[i] + '\n')
    fnew.close()

def getColourCode(spName, tempName, decayRate, decayPop):
    mlDist = calculateMaxLikDist(spName, tempName)
    #print(mlDist)
    if decayRate < 0.01:
        traceability = 1
    else:
        traceability = 1 - ((decayPop * math.exp(decayRate * mlDist)) / (1 + decayPop * (math.exp(decayRate * mlDist) - 1)))
    #print(traceability)
    if traceability <= 1 and traceability >= 0.9:
        colCode = '16711936'
    elif traceability < 0.9 and traceability >= 0.8:
        colCode = '3604736'
    elif traceability < 0.8 and traceability >= 0.7:
        colCode = '2294016'
    elif traceability < 0.7 and traceability >= 0.6:
        colCode = '983296'
    elif traceability < 0.6 and traceability >= 0.5:
        colCode = '256'
    elif traceability < 0.5:
        colCode = '65536'

    return colCode, traceability

def main(nexusTreeFile, mapFile, protId, spTree, plotFigTree, speciesMaxLikFile, speciesId, cache_dir,intended_fas_score_calc):
    global sp_max_file, hamstrMapFile, cacheDir

    matrixDict = {}
    fasFile = 'ogSeqs_' + protId + '.fasScore'
    orthFile = 'ogSeqs_' + protId + '.fa'
    sp_max_file = speciesMaxLikFile
    decayRate = 0.1
    decayPop = 0.04
    cacheDir = cache_dir

    for line in open(mapFile):
        orth = "NA"
        fas = "NA"
        currentSpecies = line.split()[-1]
        if currentSpecies not in matrixDict:
            matrixDict[currentSpecies] = []

        if currentSpecies == speciesId:
            if os.path.exists(fasFile):
                if open(fasFile).read().split('\n')[0] != "":
                    orth = open(fasFile).read().split('\n')[0].split()[1]
                fas = "1.00"
            else:
                orth = "1"

            matrixDict[currentSpecies].append(orth + '#' + fas)
        else:
            foundSpeciesFlag = False
            if os.path.exists(fasFile):
                for line2 in open(fasFile):
                    line2_species_id = line2.split()[2].split('_')[1] if '_' in line2 else line2.split()[2]
                    if line2_species_id == currentSpecies:
                        foundSpeciesFlag = True
                        orth = line2.split()[3]
                        fas = line2.split()[-1]
                        matrixDict[currentSpecies].append(orth + '#' + fas)
            elif os.path.exists(orthFile):
                orth = "0"
                for line2 in open(orthFile):
                    if '>' in line2:
                        line2_species_id = line2.split()[0].split('_')[1] if '_' in line2 else line2.split()[0][1:]
                        if line2_species_id == currentSpecies:
                            foundSpeciesFlag = True
                            orth = "1"
                            matrixDict[currentSpecies].append(orth + '#' + fas)
            if not foundSpeciesFlag:
                matrixDict[currentSpecies].append(orth + '#' + fas)

    #print('Dictionary length: ', len(matrixDict))
    try:
        hamstrMapFile = open(mapFile).read().split('\n')
        speciesTree = open(spTree).read()
        decayRate = float(open('decay_summary_%s.txt_parameter' %protId).read().split('\n')[1])
        decayPop = float(open('decay_summary_%s.txt_parameter' %protId).read().split('\n')[0])

    except IOError:
        print('ERROR: Colourizing tree encountered problem!!!')

    traceResults = open('trace_results_%s.txt' %protId, 'w')

    if os.path.exists(nexusTreeFile):
    #speciesId = nexusTreeFile.split('_')[1].split('.')[0][:5]
        for i in range(len(hamstrMapFile) - 1):
            if speciesId == hamstrMapFile[i].split('\t')[-1]:
                speciesName = hamstrMapFile[i].split('\t')[1]
                break

        colourize(speciesName, speciesId, nexusTreeFile, speciesTree, decayRate, decayPop, traceResults, matrixDict)

        try:
            subprocess.check_output('java -cp %s figtreepdf %s' %(plotFigTree, nexusTreeFile.replace('.nexus', '_edit.nexus')),shell=True)
            print("Done!")
        except subprocess.CalledProcessError as e:
            print(e.output)
            print('WARNING: No representation of traceabilities on tree possible.\nJAVA program figtreepdf not responding!!!')
    traceResults.close()

    print('Creating matrix file for PhyloProfile...')
    matrixFile = open('%s_phyloMatrix.txt' %protId, 'w')
    if intended_fas_score_calc:
        matrixFile.write('geneID\tncbiID\torthoID\tFAS\tTraceability\n')
    else:
        matrixFile.write('geneID\tncbiID\torthoID\tTraceability\n')
    for j in range(len(hamstrMapFile) - 1):
        ncbiId = "ncbi" + hamstrMapFile[j].split()[-2]
        oma_id = hamstrMapFile[j].split()[-1]
        for elements in matrixDict[oma_id]:
            if intended_fas_score_calc:
                matrixFile.write(protId + '\t' + ncbiId + '\t' + elements.split('#')[0] + '\t' + elements.split('#')[1] + '\t' + elements.split('#')[2] + '\n')
            else:
                matrixFile.write(protId + '\t' + ncbiId + '\t' + elements.split('#')[0] + '\t' + elements.split('#')[2] + '\n')
    matrixFile.close()

