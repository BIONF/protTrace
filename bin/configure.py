# Module to read in program config file and prepare variables
# to be imported in other modules

import os, sys

class setParams:

    def __init__(self, configFile):
        config = open(configFile).read().split('\n')
        for line in config:
            if line != "" and line[0] != "#":
                if line.split(':')[0] == 'species':
                    self.species = line.split(':')[1]
                if line.split(':')[0] == 'search_oma_database':
                    if line.split(':')[1] == 'YES':
                        self.search_oma_database = True
                    else:
                        self.search_oma_database = False
                if line.split(':')[0] == 'orthologs_prediction':
                    if line.split(':')[1] == 'YES':
                        self.orthologs_prediction = True
                    else:
                        self.orthologs_prediction = False
                if line.split(':')[0] == 'run_hamstr':
                    if line.split(':')[1] == 'YES':
                        self.run_hamstr = True
                    else:
                        self.run_hamstr = False
                if line.split(':')[0] == 'run_hamstrOneSeq':
                    if line.split(':')[1] == 'YES':
                        self.run_hamstrOneSeq = True
                    else:
                        self.run_hamstrOneSeq = False
                        self.run_hamstr = False
                if line.split(':')[0] == 'fas_score':
                    if line.split(':')[1] == 'YES':
                        self.fas_score = True
                    else:
                        self.fas_score = False
                if line.split(':')[0] == 'preprocessing':
                    if line.split(':')[1] == 'YES':
                        self.preprocessing = True
                    else:
                        self.preprocessing = False
                if line.split(':')[0] == 'traceability_calculation':
                    if line.split(':')[1] == 'YES':
                        self.traceability_calculation = True
                    else:
                        self.traceability_calculation = False
                if line.split(':')[0] == 'calculate_scaling_factor':
                    if line.split(':')[1] == 'YES':
                        self.calculate_scaling_factor = True
                    else:
                        self.calculate_scaling_factor = False
                if line.split(':')[0] == 'calculate_indel':
                    if line.split(':')[1] == 'YES':
                        self.calculate_indel = True
                    else:
                        self.calculate_indel = False
                if line.split(':')[0] == 'perform_msa':
                    if line.split(':')[1] == 'YES':
                        self.perform_msa = True
                    else:
                        self.perform_msa = False
                if line.split(':')[0] == 'delete_temporary_files':
                    if line.split(':')[1] == 'YES':
                        self.delete_temp = True
                    else:
                        self.delete_temp = False
                if line.split(':')[0] == 'reuse_cache':
                    if line.split(':')[1] == 'YES':
                        self.reuse_cache = True
                    else:
                        self.reuse_cache = False
                if line.split(':')[0] == 'map_traceability_tree':
                    if line.split(':')[1] == 'YES':
                        self.mapTraceabilitySpeciesTree = True
                    else:
                        self.mapTraceabilitySpeciesTree = False
                if line.split(':')[0] == 'include_paralogs':
                    if line.split(':')[1] == 'YES':
                        self.includeParalogs = True
                    else:
                        self.includeParalogs = False
                if line.split(':')[0] == 'orthologs_tree_reconstruction':
                    if line.split(':')[1] == 'YES':
                        self.phylogeneticTreeReconstruction = True
                    else:
                        self.phylogeneticTreeReconstruction = False
                if line.split(':')[0] == 'run_spartaABC':
                    if line.split(':')[1] == 'YES':
                        self.run_spartaABC = True
                    else:
                        self.run_spartaABC = False
                if line.split(':')[0] == 'dawg_instead_of_indelible':
                    if line.split(':')[1] == 'YES':
                        self.evolve_dawg = True
                    else:
                        self.evolve_dawg = False
                if line.split(':')[0] == 'aa_substitution_matrix':
                    self.aa_substitution_matrix = line.split(':')[1]
                if line.split(':')[0] == 'default_indel':
                    self.default_indel = line.split(':')[1]
                if line.split(':')[0] == 'default_indel_distribution':
                    self.default_indel_distribution = line.split(':')[1]
                if line.split(':')[0] == 'default_scaling_factor':
                    self.default_scaling_factor = line.split(':')[1]
                if line.split(':')[0] == 'simulation_runs':
                    self.simulation_runs = int(line.split(':')[1])
                if line.split(':')[0] == 'nr_of_processors':
                    self.nr_processors = int(line.split(':')[1])
                if line.split(':')[0] == 'spartaABC':
                    self.sparta = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'linsi':
                    self.msa = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'REvolver':
                    self.REvolver = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'hmmfetch':
                    self.hmmfetch = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'hmmscan':
                    self.hmmscan = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'iqtree':
                    self.iqtree = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'clustalw':
                    self.clustalw = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'blastp':
                    self.blastp = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'makeblastdb':
                    self.makeblastdb = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'Rscript':
                    self.R = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'hamstr':
                    self.hamstr = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'oneseq':
                    self.hamstrOneSeq = os.path.abspath(line.split(':')[1])
                #if line.split(':')[0] == 'degapping':
                #	self.degapping = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'species_MaxLikMatrix':
                    self.species_MaxLikMatrix = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'path_output_dir':
                    self.path_work_dir = os.path.abspath(line.split(':')[1])
                    if not os.path.exists(self.path_work_dir):
                        os.mkdir(self.path_work_dir)
                if line.split(':')[0] == 'path_cache':
                    self.path_cache = os.path.abspath(line.split(':')[1])
                    if not os.path.exists(self.path_cache):
                        os.mkdir(self.path_cache)
                if line.split(':')[0] == 'Xref_mapping_file':
                    self.hamstr_oma_tree_map = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'path_oma_seqs':
                    self.path_oma_seqs = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'path_oma_group':
                    self.path_oma_group= os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'pfam_database':
                    self.pfam_database = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'reference_species_tree':
                    self.species_tree = os.path.abspath(line.split(':')[1])
                #if line.split(':')[0] == 'species_tree_msa':
                #	self.species_tree_msa = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'simulation_tree':
                    self.simulation_tree = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'decay_script':
                    self.decay_script = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'plot_figtree':
                    self.plot_figtree = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'fas_annotations':
                    self.fas_annotations = os.path.abspath(line.split(':')[1])
                if line.split(':')[0] == 'hamstr_environment':
                    if not line.split(':')[1] == "default":
                        self.hamstr_environment = line.split(':')[1]
                    else:
                        self.hamstr_environment = ""



