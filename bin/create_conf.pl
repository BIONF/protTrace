#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use List::Util 'first';
use LWP::Simple;

# Copyright (C) 2018 INGO EBERSBERGER, ebersberger@bio.uni-frankfurt.de
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 3 of the License
# or any later version.

# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; If not, see http://www.gnu.org/licenses


#### Synopsis
##A simple perl script to check the dependencies for protTrace and for setting up the configure script

#### Default settings and dependencies ##########
my @dependencies = qw(python perl java linsi hmmscan hmmfetch makeblastdb blastp iqtree Rscript);
my @neDependencies = qw(oneseq.pl);
my $minVersion->{linsi}=6;
$minVersion->{hmmscan}=3.1;
$minVersion->{blastp}=2.7;
$minVersion->{iqtree}=1.671;
$minVersion->{Rscript} = 3;
$minVersion->{python}=2.7;
$minVersion->{java}=1.7;
$minVersion->{"oneseq.pl"}=1;
$minVersion->{figtree}=1.4;
##################
my $omaLink = 'https://omabrowser.org/All/oma-groups.txt.gz';
my $omaSeqLink = 'https://omabrowser.org/All/oma-seqs.fa.gz';
my $pfamLink = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz';
##################
my $currwd = `pwd`;
chomp $currwd;
$currwd =~ s/(.*protTrace\/).*/$1/;
#################
## some variables (values will be filled in sub userOptions)
my %prepOptions;
my %optionsValues;
my %usedFiles;
my @generalOptions;
my @preProcessing;
my @prepAdvanced;
my @scaling;
my @indel ;
my @trace;
my @path2Deps;
my @usedFileList;
my @checkList;
my @OptionList = userOptions();
my $message;
##################
my $helpmessage = "SYNOPSIS
This script checks for the essential program dependencies of protTrace and creates the config file controlling the protTrace run

USAGE: create_conf.pl -name [-update -hamstr -getOma]

OPTIONS
-name=<>	specify the name of the config file
-update		set this flag to update an existing config file provided via -name
-hamstr		set this flag if you want to make use of the HaMStR package. If not set, missing HaMStR dependencies are not traced.
-getOma		set this flag to retrieve OMA orthologs from the OMA database
-getPfam	set this flag to retrieve the Pfam-A database";
####################
## Options
my $update; ## update an existing script provided via -name
my $name; ## the name of the configure script.
my $help;
#my $install; ## flag that sets whether conda will be used to install the dependencies that are available via the conda package manager
my $full; ## set this flag to check also for the non-essential dependencies;
my $paths;
my $hamstr;
my @log = qw();
my @textarray = qw();
my $getOma;
my $getPfam;

GetOptions (
	"update" => \$update,
	"name=s" => \$name,
	"help" => \$help,
	"h"	=> \$help,
	"hamstr" => \$hamstr,
	"getOma" => \$getOma,
	"getPfam" => \$getPfam
	);
if ($help or !defined $name){
	print "\n\n\n$helpmessage\n\n";
	exit;
}
################
## create paths if not exist
if (! -e "$currwd/output") {
	my $suc = `mkdir $currwd/output`;
	if ($suc) {
		$message = "Could not create directory $currwd/output. Make sure to solve this issue before running protTrace\n";
		push @log, $message;
		print "$message\n";
	}
	else {
		$message = "Successfully created output dir $currwd/output";
		push @log, $message;
		print "$message\n";
	}
}
if (! -e "$currwd/cache") {
        my $suc = `mkdir $currwd/cache`;
        if ($suc) {
                $message = "Could not create directory $currwd/cache. Make sure to solve this issue before running protTrace\n";
                push @log, $message;
                print "$message\n";
        }
        else {
                $message = "Successfully created output dir $currwd/cache";
                push @log, $message;
                print "$message\n";
        }
}
################
($name, @textarray) = processInput($name);
if (!defined $name){
	print"No valid filename was provided. Terminating...\n";
	exit;
}
####### get the config text
my @configDat = qw();
if ($update) {
	@configDat = readConfig($name);
	populateDefault(@configDat);	
}

################
## The updating option
if (! $update){
	## first time to run the script, do the full check for dependencies
	for (my $i = 0; $i < @dependencies; $i++) {
		my $path2prog = testDep($dependencies[$i], $minVersion->{$dependencies[$i]});
		$paths->{$dependencies[$i]} = $path2prog;
		if ($paths->{$dependencies[$i]} eq '') {
			push @log, "Dependency $dependencies[$i] could not be established. Please fix this dependency before running protTest.";
			$prepOptions{$dependencies[$i]}='DEPENDENCY MISSING';
		}
		else {
			push @log, "Dependency $dependencies[$i] successfully established. Path: $paths->{$dependencies[$i]}";
			$prepOptions{$dependencies[$i]}=$paths->{$dependencies[$i]};
		}
	}
	## print output
	for (my $i = 0; $i < @dependencies; $i++) {
			print "Dependency successfully established: ";
			print "$paths->{$dependencies[$i]}\n";
	}
}
if (1){
	my @updateKeys = qw();
	my @select;
	my $entry;
	my @toCheck = getCheckList();
	if (!defined $toCheck[0]){
		print "nothing to update\n";
		push @log, "nothing to update";
	}
	else {
		my $count = 0;
		for (my $i = 0; $i < @toCheck; $i++){			
			print "$checkList[$toCheck[$i]]\n";
			$entry = $OptionList[$toCheck[$i]];
			for (my $j = 0; $j < @$entry; $j++){
			print "\t[$count] - @$entry[$j]\n";
			$select[$count] = @$entry[$j];
			$count ++;
			}
			print "\n";
			@updateKeys = (@updateKeys, userInput("numeric", "Select the entries for updating or type q to quit: "));
			if (!defined $updateKeys[0]){
				print "nothing to update\n";
				push @log, "nothing to update";
				printLog();
			}
		}
		if (@updateKeys){
			for (my $i = 0; $i < @updateKeys; $i++){
				my $message = "$i - $select[$updateKeys[$i]]\n\tCurrent value: $prepOptions{$select[$updateKeys[$i]]}";
				print "$message\n";
				my $message2 = "\tEnter new value or type k to keep old [$optionsValues{$select[$updateKeys[$i]]}]: ";
				push @log, $message2;
				my $newVal=userInput("alpha", $message2);
				if (defined $newVal and $newVal ne '' and $newVal ne 'k'){
					if ($newVal =~ /\//){
						## The value is a path and we check for its existence
						if (! -e $newVal){
							push @log, "You entered a path that does not exist. The configure script will continue but you need to fix this before running protTrace";
							print "You entered a path that does not exist. The configure script will continue but you need to fix this before running protTrace\n\n";
						}
						else {
							push @log, "You provided $newVal, path exists";
							print "You provided $newVal, path exists\n\n";
						}
					}
					$prepOptions{$select[$updateKeys[$i]]} = $newVal;
				}
				elsif (defined $newVal and $newVal eq 'k') {
					my $message = "Nothing to update, keeping old";
					print "$message\n";
					push @log, $message;
				}
				elsif (defined $newVal and $newVal eq ''){
					$newVal = $optionsValues{$select[$updateKeys[$i]]};
					$newVal =~ s/.*\(Default (.*)\)/$1/; 
					my $message = "Using default value: $newVal\n";
					push @log, $message;
					print "$message\n";
					$prepOptions{$select[$updateKeys[$i]]} = $newVal; 
				}
					
				else {
					my $message = "You chose to quit";
					print $message . "\n";
					push @log, $message;
					printLog();
					exit;
				}
			}
		}
	}
}
checkPaths();
############## get the OMA orthologs
if ($getOma){
	my $message;
	my $omaret = retrieveData($omaLink, $prepOptions{path_oma_group});
	my $omaret2= retrieveData($omaSeqLink, $prepOptions{path_oma_seqs});
	if ($omaret && $omaret2) {
		$message = "Oma files are in place under $prepOptions{path_oma_group} and $prepOptions{path_oma_seqs}";
		push @log, $message;
		print "$message\n";
                print "Reformatting fasta\n";
                my $return = reformatFasta($prepOptions{path_oma_seqs});
		if ($return) {
			$message = "Something went wrong during reformatting of the oma sequences. Please check\n";
			push @log, $message;
			print "$message\n";
		}

	}
	else {
		$message = "Problems during retrieval of OMA data. Files under $prepOptions{path_oma_group} and $prepOptions{path_oma_seqs} are not in place. Please solve before running protTrace";
		push @log, $message;
		print "$message\n";
	}
}
############# get the Pfam data
if ($getPfam) {
	my $message;
	my $pfamret = retrieveData($pfamLink, $prepOptions{pfam_database});
	if ($pfamret) {
		$message = "running hmmpress on $prepOptions{pfam_database}\n";
		push @log, $message;
		my $pressRet = `hmmpress $prepOptions{pfam_database}`;
		if ($pressRet){
			$message = "Pfam HMMs are in place under $prepOptions{pfam_database}";
		}
		else {
			$message = "Pfam-A downloaded to $prepOptions{pfam_database}, but hmmpress failed. Please solve before running protTrace";
		}
                push @log, $message;
                print "$message\n";
        }
        else {  
                $message = "Problems during retrieval of Pfam data. Files under $prepOptions{pfam_database}. Please solve before running protTrace";
                push @log, $message;
                print "$message\n";
        }
}

printConfig($name);
printLog();
exit;

############### Start Sub
sub testDep {
	my ($prog, $minversion) = @_;
	my $jump;
	my $exist;
	my $path2prog;
	my $issue;
	print "Testing for $prog\n";
	while (! $exist and ! $jump) {
		$path2prog = `which $prog`;
		chomp $path2prog;
		if ($path2prog =~ /\//) {
			$exist = 1;
			## version check only for the critical programs
			my $version;
			my $versionshort;
			if ($prog eq 'iqtree' or $prog eq 'Rscript' or $prog eq 'python') {
				$version = `$prog --version 2>&1 |head -n 1`;
				chomp $version;
				$version =~ s/.* ([0-9.]{1,}).*/$1/;
				$versionshort = $version;
				$versionshort =~ s/([0-9]{1,}).([0-9.]{1,})/$1_$2/g;
				$versionshort =~ s/\.//g;
				$versionshort =~ s/_/./;
				if ($versionshort < $minversion){
					## problem
					my $message = "Issue: I found $prog, but the version number $versionshort is below the minimal requirement $minversion. You can continue with the configure script, but you will have to upgrade $prog manually";
					push @log, $message;
					print "$message\n";
					$issue = 1;
				}
				if ($prog eq 'python' and $versionshort > 3) {
					my $message = "Fatal: I found python version $versionshort but protTrace requires python 2.7. Please install python 2.7 and restart the script";
					print "$message\n";
					printLog();
					exit;
				}
			} 
				  
		}
		else {
			$issue = 1;
		}
		if ($issue){
			my $message = "There is an issue with $prog. Please select [p(path)|q(uit)|s(kip)]: ";
			print $message;
			push @log, $message;
			my $q=<>;
			if ($q =~ /^q/i) {
				print "you chose to quit. Exiting...\n";
				printLog();
				exit;
			}
			elsif ($q =~ /^p/i) {
				while (! -e $prog){
					my $message = "I could not find $prog! Please provide the absolute path to $prog or type q to quit: ";
					print "$message";
					push @log, $message;
					$path2prog=<STDIN>;
					chomp $path2prog;
					if ($path2prog eq 'q'){
						my $message = "You chose to quit. Exiting...\n";
						push @log, $message;
						printLog();
						exit;
					}
					$path2prog =~ s/\/$//;
					if ($path2prog !~ /$prog$/) {
						$path2prog = $path2prog . '/' .  $prog;
						$prog = $path2prog;
					}
				}
				$exist = 1;				
			}
			elsif ($q =~ /^s/i) {
				$jump = 1;
				my $message = "You chose to skip the test for $prog. The setup script will continue, but protTrace will not be in a runable state unless you install the missing dependency of $prog manually";
				push @log, $message;
				print "$message\n";
				$path2prog = '';
			}
		}
	}
	return ($path2prog);
}	
				 
###################
## get user input
sub userInput {
	my ($type, $message) = @_;
	if (defined $message) {
		print "$message";
	} 
	my $answer;
	$answer = <STDIN>;
	chomp $answer;
	if ($answer =~ /^q/i){
		return();
	}
	if ($type eq 'numeric'){
		while ($answer =~ /[^0-9,-]/ and $answer !~ /^q/i){
			print "please enter only a (comma separated) list of digits or a range: ";
			$answer = <STDIN>;
			$answer =~ s/ //g;
			chomp $answer;
		}
		if ($answer =~ /^q/i){
			return();
		}
		else {
			my @list = split ',', $answer;
			my @listret = qw();
			for (my $i = 0; $i < @list; $i++){
				if ($list[$i] !~ /-/){
					push @listret, $list[$i];
				}
				elsif ($list[$i] =~ /([0-9]+)-([0-9]+)/){
					for (my $j = $1; $j <= $2; $j++){
						push @listret, $j;
					}
				}
			}
			return(@listret); 
		}
	}
	elsif ($type eq 'alpha'){
		return($answer); 
	}
}
################################
## testing for the user input
sub processInput {
	my $nameSub = shift;
	while (! defined $nameSub) {
		print "Please enter a name for the config file or type q for quitting:";
		$nameSub=<STDIN>;
		chomp $nameSub;
		if ($nameSub =~ /^q$/){
			return();
		}
	}
	if (defined $update && ! (-e $nameSub)) {
		my $message = "You specified the option -update but the file $name seems to not exist in $currwd. You can provide the path or enter (q) to cancel: "; 
		push @log, $message;
		print "$message\n";
		my $stop = 0;
		while ($stop == 0 and ! (-e $nameSub)){
			my $answer = <STDIN>;
			chomp $answer;
			push @log, "\t$answer";
			if (!defined $answer or $answer eq '' or $answer =~ /^q/i){
				$stop = 1;
				$nameSub = undef;
			}
			elsif ($answer !~ /$nameSub/){
				$nameSub = $answer . '/' . $nameSub;
			}
			else {
				$nameSub = $answer;
			}
		}
	}
	return($nameSub);
}
##################
sub populateDefault {
	for (keys %prepOptions){
		$prepOptions{$_} = readValue($_,@configDat);
	}
	return();
}

##################
sub readValue {
	my ($part, @list) = @_;
	my $value = first { /^$part/ } @configDat;
	$value =~ s/^$part://;
	chomp $value;
	return($value);
}
##################
sub getCheckList {
	@checkList = ("General Options", "Preprocessing Settings", "Advanced Preprocessing Settings", "Scaling factors",
	"Indel parameter", "Traceability calculation", "Program paths", "Paths to files");	
	print "\nPlease select from the list the groups you wish to update or enter q to quit. To use defaults, just hit enter:\n";
	for (my $i = 0; $i < @checkList; $i++){
		print "[$i] - $checkList[$i]\n";
	}
	my @list = userInput("numeric", "If you select more than one group, separate them by a comma or give a range. Leave blank to initialize with default values: ");
	return(@list); 
}
###################
sub userOptions {
	%prepOptions = (species => 'YEAST',
					  nr_of_processors => 1,
					  delete_temporary_files => 'YES',
					  reuse_cache => 'YES',
					  preprocessing => 'YES',
					  orthologs_prediction => 'YES',
					  search_oma_database => 'YES',
					  run_hamstr => 'NO',
					  run_hamstrOneSeq => 'NO',
					  include_paralogs => 'YES',
					  fas_score => 'NO',
					  orthologs_tree_reconstruction => 'YES',
					  calculate_scaling_factor => 'YES',
					  default_scaling_factor => 1,
					  perform_msa => 'YES',
					  calculate_indel => 'YES',
					  default_indel => 0.08,
					  default_indel_distribution => 0.25,
					  traceability_calculation => 'YES',
					  aa_substitution_matrix => 'WAG',
					  simulation_runs => 100,
					  path_output_dir => "$currwd/output",
					  path_cache => "$currwd/cache",
					  map_traceability_tree => 'YES',
					  REvolver => "$currwd/used_files/REvolver.jar",
					simulation_tree => "$currwd/used_files/stepWiseTree.nw",
					decay_script => "$currwd/used_files/r_nonlinear_leastsquare.R",
					plot_figtree => "$currwd/used_files/plotPdf.jar",
					Xref_mapping_file => "$currwd/used_files/speciesTreeMapping.txt",
					reference_species_tree => "$currwd/used_files/speciesTree.nw",
					species_MaxLikMatrix => "$currwd/used_files/speciesLikelihoodMatrix.txt",
					path_oma_seqs => "$currwd/used_files/oma-seqs.fa",
					path_oma_group => "$currwd/used_files/oma-groups.txt",
					pfam_database => "$currwd/used_files/Pfam-A.hmm",
					fas_annotations => 'PROVIDEPATH2HaMStR/weight_dir',
					hamstr_environment => 'default',
					iqtree => "$currwd/used_files/iqtree",
					linsi => '/share/applications/bin/linsi',
					hmmfetch => '/share/applications/bin/hmmfetch',
					hmmscan => '/share/applications/bin/hmmscan',
					blastp => '/share/applications/bin/blastp',
					makeblastdb => '/share/applications/bin/makeblastdb',
					Rscript => '/usr/bin/Rscript',
					hamstr => '',
					oneseq =>  '');		
					
	%optionsValues = (species => 'OMA 5 letter Species Identifier: (Default YEAST)',
					  nr_of_processors => 'Integer (Default 1)',
					  delete_temporary_files => 'YES|NO (Default NO)',
					  reuse_cache => 'YES|NO (Default YES)',
					  preprocessing => 'YES|NO (Default YES)',
					  orthologs_prediction => 'YES|NO (Default YES)',
					  search_oma_database => 'YES|NO Default YES)',
					  run_hamstr => 'YES|NO (Default NO)',
					  run_hamstrOneSeq => 'YES|NO (Default NO)',
					  include_paralogs => 'YES|NO (Default NO)',
					  fas_score => 'YES|NO (Default NO)',
					  orthologs_tree_reconstruction => 'YES|NO (Default YES)',
					  calculate_scaling_factor => 'YES|NO (Default YES)',
					  default_scaling_factor => 'float (Default 1.57)',
					  perform_msa => 'YES|NO (Default YES)',
					  calculate_indel => 'YES|NO (Default YES)',
					  default_indel => 'float (Default 0.08)',
					  default_indel_distribution => 'float (Default 0.25)',
					  traceability_calculation => 'YES|NO (Default YES)',
					  aa_substitution_matrix => 'JTT|WAG|LG|Blosum62|mtMAM|mtREV|mtART (Default WAG)',
					  simulation_runs => 'Integer (Default 100)',
					  path_output_dir => "Path to output directory (Default $currwd/output)",
					  path_cache => "Path to cache directory (Default $currwd/cache)",
					  map_traceability_tree => 'YES|NO (Default NO)',
					  REvolver => "Path to REvolver (Default $currwd/used_files/REvolver.jar)",
					simulation_tree => "Path to simulation tree (Default $currwd/used_files/stepWiseTree.nw)",
					decay_script => "Path to decay script (Default $currwd/used_files/r_nonlinear_leastsquare.)",
					plot_figtree => "Path to tree plotter (Default $currwd/used_files/plotPdf.jar)",
					Xref_mapping_file => "Path to mapping file (Default $currwd/used_files/speciesTreeMapping.txt)",
					reference_species_tree => "Path to species tree (Default $currwd/used_files/speciesTree.nw)",
					species_MaxLikMatrix => "Path to distance matrix (Default $currwd/used_files/speciesLikelihoodMatrix.txt)",
					path_oma_seqs => "Path to OMA sequences (Default $currwd/used_files/oma-seqs.fa)",
					path_oma_group => "Path to OMA groups (Default $currwd/used_files/oma-groups.txt)",
					pfam_database => "Path to Pfam (Default $currwd/used_files/Pfam-A.hmm)",
					fas_annotations => 'Path to HaMStR weight_dir (Default NULL)',
					hamstr_environment => 'Path to HaMStR directory (Default NULL)',
					iqtree => 'Path to iqtree (Default NULL)',
					linsi => 'Path to linsi (Default NULL)',
					hmmfetch => 'Path to hmmfetch (Default NULL)',
					hmmscan => 'Path to hmmscan (Default NULL)',
					blastp => 'Path to blastp (Default NULL)',
					makeblastdb => 'Path to makeblastdb (Default NULL)',
					Rscript => 'Path to Rscript (Default NULL)',
					hamstr => 'Path to hamstr.pl (Default NULL)',
					oneseq =>  'Path to oneseq.pl (Default NULL)');						
	@generalOptions = qw(species nr_of_processors delete_temporary_files reuse_cache map_traceability_tree);
	@preProcessing = qw(preprocessing orthologs_prediction search_oma_database);
	@prepAdvanced = qw(run_hamstr run_hamstrOneSeq include_paralogs fas_score);
	@scaling = qw(calculate_scaling_factor default_scaling_factor);
	@indel = qw(perform_msa calculate_indel default_indel default_indel_distribution);
	@trace = qw(traceability_calculation aa_substitution_matrix simulation_runs);
	@path2Deps = qw(iqtree linsi hmmfetch hmmscan blastp makeblastdb Rscript hamstr oneseq);
	@usedFileList = qw(REvolver simulation_tree decay_script plot_figtree Xref_mapping_file reference_species_tree
	species_MaxLikMatrix path_oma_seqs path_oma_group pfam_database fas_annotations hamstr_environment path_output_dir path_cache); 

	my @options;
	$options[0] = \@generalOptions;
	$options[1] = \@preProcessing;
	$options[2] = \@prepAdvanced;
	$options[3] = \@scaling;
	$options[4] = \@indel;
	$options[5] = \@trace;
	$options[6] = \@path2Deps;
	$options[7] = \@usedFileList;
	return(@options);	
}

##########################
sub checkPaths {
	for (my $i = 0; $i < @usedFileList; $i++){
		if (! -e $prepOptions{$usedFileList[$i]}){
			if (($usedFileList[$i] =~ /hamstr/i or $usedFileList[$i] =~ /fas_annotations/) and ! $hamstr){
				next;
			}
			else {
				my $message = "The files under $prepOptions{$usedFileList[$i]} does not exist. The configure script will run to the end, but you will need to fix this before running protTrace";
				push @log, $message;
				print "$message\n";
			}
		}
	}
}
##########################
sub readConfig {
	my $nameSub2 = shift;
	open (IN, "$nameSub2") or die "could not open $nameSub2 for reading\n";
	my @configText = <IN>;
	close IN or die "could not close Filenhandle for $nameSub2 after reading\n";
	return(@configText);
}
#########################
sub printConfig {
	my $name = shift;
	open (OUT, ">$name") or die "could not open $name for writing output\n";
	print OUT "#####   protTrace Run Configuration     #####\n";
	print OUT "species:$prepOptions{species}\n";
	print OUT "nr_of_processors:$prepOptions{nr_of_processors}\n";
	print OUT "delete_temporary_files:$prepOptions{delete_temporary_files}\n";
	print OUT "reuse_cache:$prepOptions{reuse_cache}\n";
	print OUT "###\n";
	print OUT "#####   Preprocessing Steps Configuration (Only valid when preprocessing is YES)        #####\n";
	print OUT "preprocessing:$prepOptions{preprocessing}\n";
	print OUT "###     Orthologs Search Steps Configuration (Only valid when orthologs_prediction is YES)      ###\n";
	print OUT "orthologs_prediction:$prepOptions{orthologs_prediction}\n";
	print OUT "search_oma_database:$prepOptions{search_oma_database}\n";
	print OUT "run_hamstr:$prepOptions{run_hamstr}\n";
	print OUT "run_hamstrOneSeq:$prepOptions{run_hamstrOneSeq}\n";
	print OUT "include_paralogs:$prepOptions{include_paralogs}\n";
	print OUT "###     FAS Score Calculations  ###\n";
	print OUT "fas_score:$prepOptions{fas_score}\n";
	print OUT "###     Orthologs Tree Reconstruction   ###\n";
	print OUT "orthologs_tree_reconstruction:$prepOptions{orthologs_tree_reconstruction}\n";
	print OUT "###     Scaling Factor Calculation      ###\n";
	print OUT "calculate_scaling_factor:$prepOptions{calculate_scaling_factor}\n";
	print OUT "default_scaling_factor:$prepOptions{default_scaling_factor}\n";
	print OUT "###     Multiple Sequence Alignment     ###\n";
	print OUT "perform_msa:$prepOptions{perform_msa}\n";
	print OUT "###     Calculate Indel rates   ###\n";
	print OUT "calculate_indel:$prepOptions{calculate_indel}\n";
	print OUT "default_indel:$prepOptions{default_indel}\n";
	print OUT "default_indel_distribution:$prepOptions{default_indel_distribution}\n";
	print OUT "#####\n";
	print OUT "###\n";
	print OUT "#####   Traceability Calculation Configuration  #####\n";
	print OUT "traceability_calculation:$prepOptions{traceability_calculation}\n";
	print OUT "aa_substitution_matrix:$prepOptions{aa_substitution_matrix}\n";
	print OUT "simulation_runs:$prepOptions{simulation_runs}\n";
	print OUT "#####\n";
	print OUT "###\n";
	print OUT "#####   Writing traceability results to file and on reference species tree      #####\n";
	print OUT "map_traceability_tree:$prepOptions{map_traceability_tree}\n";
	print OUT "#####\n";
	print OUT "###\n";
	print OUT "#####   Configuring paths for protTrace dependencies    #####\n";
	print OUT "iqtree:$prepOptions{iqtree}\n";
	print OUT "linsi:$prepOptions{linsi}\n";
	print OUT "hmmfetch:$prepOptions{hmmfetch}\n";
	print OUT "hmmscan:$prepOptions{hmmscan}\n";
	print OUT "blastp:$prepOptions{blastp}\n";
	print OUT "makeblastdb:$prepOptions{makeblastdb}\n";
	print OUT "Rscript:$prepOptions{Rscript}\n";
	print OUT "hamstr:$prepOptions{hamstr}\n";
	print OUT "oneseq:$prepOptions{oneseq}\n";
	print OUT "\n###\n";
	print OUT "#####   Used files / directories #####\n";
	print OUT "REvolver:$prepOptions{REvolver}\n";
	print OUT "simulation_tree:$prepOptions{simulation_tree}\n";
	print OUT "decay_script:$prepOptions{decay_script}\n";
	print OUT "plot_figtree:$prepOptions{plot_figtree}\n";
	print OUT "Xref_mapping_file:$prepOptions{Xref_mapping_file}\n";
	print OUT "reference_species_tree:$prepOptions{reference_species_tree}\n";
	print OUT "species_MaxLikMatrix:$prepOptions{species_MaxLikMatrix}\n";
	print OUT "path_oma_seqs:$prepOptions{path_oma_seqs}\n";
	print OUT "path_oma_group:$prepOptions{path_oma_group}\n";
	print OUT "pfam_database:$prepOptions{pfam_database}\n";
	print OUT "fas_annotations:$prepOptions{fas_annotations}\n";
	print OUT "hamstr_environment:$prepOptions{hamstr_environment}\n";
	print OUT "\n####################\n";
	print OUT "###     NOTE:- If no hamstr_environment is used or default HaMStR directories are used, then simply write hamstr_environment:default
####################\n";
	print OUT "###\n";
	print OUT "#####   Path Configuration (Where outputs will be saved) - Change this only if you want output and cache results at a different place then default. #####\n";
	print OUT "path_output_dir:$prepOptions{path_output_dir}\n";
	print OUT "path_cache:$prepOptions{path_cache}\n";
	close OUT or die "could not close OUT after writing config file $name\n";
}
##############
sub printLog{
	open OUT, ">config.log" or die "could not open logfile config.log for writing\n";
	print OUT join "\n", @log;
	print OUT "\n";
	close OUT or die "could not properly close filehandle of log file\n";
}
##############
sub retrieveData {
	my ($URL, $destination) = @_;
	my $message = "retrieving data from $URL and copying to $destination";
	push @log, $message;
	print $message . "\n";
	my $answer = getAnswerFile($destination);
	if ($answer =~ /^y/i){
		getstore($URL, "$destination.gz");
		push @log, "Downloading $destination";
		print "unzipping $destination.gz\n";
		`gunzip "$destination.gz"`;
	}
	else {
		$message = "$destination already exists, and you chose to keep it";
		push @log, $message;
		print "$message\n";
	}
	if (-e $destination) {
		return(1);
	}
	else {
		return (0);
	}
}
###################
sub getAnswerFile {
	my $filename = shift;
	my $answer = 'y';
	my $message = '';
	if (-e $filename){
		$message = "file $filename already exists! Do you want to overwrite?[y|n]: ";
		push @log, $message;
		print $message;
		$answer = '';
		while ($answer !~ /^[yn]/i){
			$answer = <>;
			if ($answer !~ /^[yn]/i){
				print "please answer with either y(es) or n(o): ";
			}
			else {
				$log[(scalar@log)-1] .= $answer;
			}
		}
	}
	return($answer);
}
###################
sub reformatFasta {
	my $fileName = shift;
	open (IN, "$fileName") or die "could not open file $fileName for reformatting\n";
	my $head;
	my $seq = '';
	my @out;
	my $count = 0;
	## remove an tmp file if it exists
	if (-e "$fileName.tmp"){
		print "Found an existing tmp file $fileName.tmp. Removing it...\n";
		`rm -f $fileName.tmp`;
	}
	while (<IN>){
		chomp $_;
		if ($_ =~ /^>/) {
			if (length($seq) > 0){
				push @out, $head;
				$head = '';
				push @out, $seq;
				$seq = '';
			}
			$head = $_;
			$head =~ s/ //g;
		}
		else {
			$seq .=$_;
		}
		if (scalar(@out) == 100000) {
			$count += scalar(@out)/2;
			open (OUT, ">>$fileName.tmp") or die "could not open tmp file for writing\n";
			print OUT join "\n", @out;
			print "$count sequences processed\n";
			@out = qw();
			close OUT or die "Could not close tmp file after writing\n";
		}
	}
	open (OUT, ">>$fileName.tmp") or die "could not open tmp file for writing\n";
        print OUT join "\n", @out;
	$count += scalar(@out)/2;
	my $message = "A total of $count sequences have been reformatted";
	print "$message\n";
	push @log, $message;
	@out = qw();
        close OUT or die "Could not close tmp file after writing\n";
	my $succ = `mv $fileName.tmp $fileName`;
	return($succ);
}
