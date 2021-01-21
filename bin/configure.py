# Module to read in program config file and prepare variables
# to be imported in other modules

import os, sys

class setParams:
    """ Stores the configuration for ProtTrace. """

    def parse_option(self, line_name_value_pair, name, parse_func):
        """ Checks the first field of the name value pair 
            for the correct name and parses the second 
            field properly. """
        if name in line_name_value_pair:
            return parse_func(line_name_value_pair, name)
        else:
            return None

    def string_option(self, line_name_value_pair, name):
        """ Parses the value of the name value pair
            to the searched name. """
        return line_name_value_pair[name]

    def integer_option(self, line_name_value_pair, name):
        """ Parses the value of the name value pair
            to the searched name into an integer. """
        return int(line_name_value_pair[name])

    def float_option(self, line_name_value_pair, name):
        """ Parses the value of the name value pair
            to the searched name into a float. """
        return float(line_name_value_pair[name])

    def boolean_option(self, line_name_value_pair, name):
        """ Checks the name value pair for YES and NO and 
            returns True and False, depending on the 
            result. """
        if line_name_value_pair[name] in set(['YES', 'yes', 'True', 'true']):
            return True
        else:
            return False

    def path_option(self, line_name_value_pair, name):
        """ Returns the absolute path to the path given
            in the name value pair. """
        return os.path.abspath(line_name_value_pair[name])

    def create_path_option(self, line_name_value_pair, name):
        """ Creates the path specified in the name value pair 
            if it does not exist yet. """
        created_path = self.path_option(line_name_value_pair, name)
        if not os.path.exists(created_path):
            os.mkdir(created_path)
        return created_path
        
    def read_config_file(self, config_file):
        """ Reads the config file in, splits it into lines,
            ignores empty and commented lines and splits each
            line into a dict of key value pairs by ':' """
        with open(config_file, 'r') as config:
            return {splitted_line[0] : splitted_line[1] for splitted_line in [line.split(':') for line in config.read().split('\n') if not line == '' and line[0] != '#']}

    def __init__(self, config_file):
        
        config_key_value_pairs = self.read_config_file(config_file)
        
        # This option sets the query species
        self.species = self.parse_option(config_key_value_pairs, 'species', self.string_option)
        # Yes-No settings
        self.search_oma_database = self.parse_option(config_key_value_pairs, 'search_oma_database', self.boolean_option)
        self.orthologs_prediction = self.parse_option(config_key_value_pairs, 'orthologs_prediction', self.boolean_option)
        self.run_hamstr = self.parse_option(config_key_value_pairs, 'run_hamstr', self.boolean_option)
        self.run_hamstrOneSeq = self.parse_option(config_key_value_pairs, 'run_hamstrOneSeq', self.boolean_option)
        if not self.run_hamstrOneSeq:
            self.run_hamstr = False
        self.fas_score = self.parse_option(config_key_value_pairs, 'fas_score', self.boolean_option)
        self.preprocessing = self.parse_option(config_key_value_pairs, 'preprocessing', self.boolean_option)
        self.traceability_calculation = self.parse_option(config_key_value_pairs, 'traceability_calculation', self.boolean_option)
        self.calculate_scaling_factor = self.parse_option(config_key_value_pairs, 'calculate_scaling_factor', self.boolean_option)
        self.calculate_indel = self.parse_option(config_key_value_pairs, 'calculate_indel', self.boolean_option)
        self.perform_msa = self.parse_option(config_key_value_pairs, 'perform_msa', self.boolean_option)
        self.delete_temp = self.parse_option(config_key_value_pairs, 'delete_temporary_files', self.boolean_option)
        self.reuse_cache = self.parse_option(config_key_value_pairs, 'reuse_cache', self.boolean_option)
        self.mapTraceabilitySpeciesTree = self.parse_option(config_key_value_pairs, 'map_traceability_tree', self.boolean_option)
        self.includeParalogs = self.parse_option(config_key_value_pairs, 'include_paralogs', self.boolean_option)
        self.phylogeneticTreeReconstruction = self.parse_option(config_key_value_pairs, 'orthologs_tree_reconstruction', self.boolean_option)
        self.run_spartaABC = self.parse_option(config_key_value_pairs, 'run_spartaABC', self.boolean_option)
        self.evolve_dawg = self.parse_option(config_key_value_pairs, 'dawg_instead_of_indelible', self.boolean_option)

        # Options where numbers or strings can be specified
        self.aa_substitution_matrix = self.parse_option(config_key_value_pairs, 'aa_substitution_matrix', self.string_option)
        self.simulation_runs = self.parse_option(config_key_value_pairs, 'simulation_runs', self.integer_option)
        self.nr_processors = self.parse_option(config_key_value_pairs, 'nr_of_processors', self.integer_option)
        
        # Default evolution simulation parameters
        self.default_indel = self.parse_option(config_key_value_pairs, 'default_indel', self.float_option)
        self.default_indel_distribution = self.parse_option(config_key_value_pairs, 'default_indel_distribution', self.float_option)
        self.default_scaling_factor = self.parse_option(config_key_value_pairs, 'default_scaling_factor', self.float_option)

        # Software paths
        self.sparta = self.parse_option(config_key_value_pairs, 'spartaABC', self.path_option)
        self.msa = self.parse_option(config_key_value_pairs, 'linsi', self.path_option)
        self.REvolver = self.parse_option(config_key_value_pairs, 'REvolver', self.path_option)
        self.hmmfetch = self.parse_option(config_key_value_pairs, 'hmmfetch', self.path_option)
        self.hmmscan = self.parse_option(config_key_value_pairs, 'hmmscan', self.path_option)
        self.iqtree = self.parse_option(config_key_value_pairs, 'iqtree', self.path_option)
        self.treepuzzle = self.parse_option(config_key_value_pairs, 'treepuzzle', self.path_option)
        self.clustalw = self.parse_option(config_key_value_pairs, 'clustalw', self.path_option)
        self.blastp = self.parse_option(config_key_value_pairs, 'blastp', self.path_option)
        self.makeblastdb = self.parse_option(config_key_value_pairs, 'makeblastdb', self.path_option)
        self.R = self.parse_option(config_key_value_pairs, 'Rscript', self.path_option)

        # HaMStR
        self.hamstr = self.parse_option(config_key_value_pairs, 'hamstr', self.path_option)
        self.hamstrOneSeq = self.parse_option(config_key_value_pairs, 'oneseq', self.path_option)
        self.hamstr_environment = self.parse_option(config_key_value_pairs, 'hamstr_environment', self.string_option)
        if self.hamstr_environment == 'default':
            hamstr_environment = ''

        # Set output and cache directory paths or create default ones
        # Each protein creates its own subdirectory within this output directory
        self.path_work_dir = self.parse_option(config_key_value_pairs, 'path_output_dir', self.create_path_option)
        # The cache directory contains global intermediary files,
        # such as calculated species distances
        self.path_cache = self.parse_option(config_key_value_pairs, 'path_cache', self.create_path_option)
        # Pairwise species distances are calculated within this directory
        self.path_distance_work_dir = self.parse_option(config_key_value_pairs, 'path_distance_work_dir', self.create_path_option)

        # Input file paths
        self.species_MaxLikMatrix = self.parse_option(config_key_value_pairs, 'species_MaxLikMatrix', self.path_option)
        self.hamstr_oma_tree_map = self.parse_option(config_key_value_pairs, 'Xref_mapping_file', self.path_option)
        self.species_MaxLikMatrix = self.parse_option(config_key_value_pairs, 'species_MaxLikMatrix', self.path_option)
        self.path_oma_seqs = self.parse_option(config_key_value_pairs, 'path_oma_seqs', self.path_option)
        self.path_oma_pairs = self.parse_option(config_key_value_pairs, 'path_oma_pairs', self.path_option)
        self.path_oma_group = self.parse_option(config_key_value_pairs, 'path_oma_group', self.path_option)
        self.pfam_database = self.parse_option(config_key_value_pairs, 'pfam_database', self.path_option)
        self.species_tree = self.parse_option(config_key_value_pairs, 'reference_species_tree', self.path_option)
        self.simulation_tree = self.parse_option(config_key_value_pairs, 'simulation_tree', self.path_option)
        self.concat_alignments_script = self.parse_option(config_key_value_pairs, 'concat_alignments', self.path_option)
        self.decay_script = self.parse_option(config_key_value_pairs, 'decay_script', self.path_option)
        self.plot_figtree = self.parse_option(config_key_value_pairs, 'plot_figtree', self.path_option)
        self.fas_annotations = self.parse_option(config_key_value_pairs, 'fas_annotations', self.path_option)
