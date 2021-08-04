#!/bin/bash
#===============================================================================
#
#          FILE: setup.sh
#
#         USAGE: ./setup.sh
#
#   DESCRIPTION: Setup script for non-python dependencies of ProtTrace.
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dominik Perisa (DP), dominik.perisa@bio.gmail.com
#  ORGANIZATION: Goethe University Frankfurt
#       CREATED: 02/10/2021 02:55:46 PM
#      REVISION:  ---
#===============================================================================


# Since the oma pairwise ortholog splits only extract species which are specified in the species mapping of ProtTrace, it can make sense to run this script again after the user changed the list.

# You need to remove or rename the original data for a rerun. This script will not replace existing data!

### Startup ###

if [ ! -d "./prottrace/" ];then
	echo "Error: Please move to the ProtTrace main directory!"
	exit
fi

u_files="$PWD/used_files"

exit_split_pairs()
{
	if [ ! -z "$split_pairs_dir" ]
	then
		if [ -d "$split_pairs_dir" ]
		then
			echo "Removing unfinished directory ${split_pairs_dir}"
			rm -r "$split_pairs_dir"
		fi
	fi
}	# ----------  end of function exit_split_pairs  ----------


exit_split_seqs()
{
	if [ ! -z "$split_seqs_dir" ]
	then
		if [ -d "$split_seqs_dir" ]
		then
			echo "Removing unfinished directory ${split_seqs_dir}"
			rm -r "$split_seqs_dir"
		fi
	fi

	exit	
}	# ----------  end of function exit_split_seqs  ----------

split_file ()
{
	source_file=$1
	split_dir=$2

	if [ ! -d "$split_dir" ]
	then
		if [ -f "$source_file" ]
		then
			echo "Created the directory to split files into."
			mkdir "$split_dir"
			return
		fi
	fi

	false 
	return
}	# ----------  end of function split_file  ----------

split_pairs ()
{
	# The pairwise orthologs file of the OMA database is very large. Therefore, we try to reduce the space and time needed by filtering the rows for species specified in the species_mapping file.
	echo "Splitting orthologous pair files. Interrupting the process with SIGTERM once will cleanup leftover files."
	if split_file $pairs_file $split_pairs_dir  
	then
		# Putting "ALL" into the species field of the configuration file enables the all-versus-all computation of species distances. This is used internally by ProtTrace, but here one needs to insert the ALL keyword beforehand to remove the filtering for a particular species.
		focal_species=$(awk 'BEGIN{FS=":";OFS=" "} $1=="species"{print $2;exit}' prog.config)

		awk -v focal=$focal_species -v pairwise_dir="$split_pairs_dir" 'BEGIN{FS=OFS="\t"} FNR==NR&&$4!=focal{arr[$4]=1;next} substr($0,1,1)=="#"{next} {s1=substr($1,1,5); s2=substr($2,1,5)} (s1==focal||s2==focal||focal=="ALL")&&(arr[s1]||arr[s2]){print $1,$2 > pairwise_dir"/"s1"_"s2".txt"}' "$u_files/species_mapping.txt" "$pairs_file"
	fi

}	# ----------  end of function split_pairs  ----------

split_seqs ()
{

	# The sequence file is separated, based on species IDs in each sequence header. By temporarily removing the header line, we can remove newlines from the sequence. Therefore, we can save the preprocessing step.
	echo "Splitting sequence files. Interrupting the process with SIGTERM once will cleanup leftover files."
	if split_file $seqs_file $split_seqs_dir
	then
		# Extract the species specified in the species_mapping file as a filter. Write them into a FASTA-like format to make them instantly parsable by the splitting function.
		awk 'BEGIN{FS="\t";OFS="\n";RS="\n";ORS=">"} {print $4}' "$u_files/species_mapping.txt" > "$u_files/sequence_extractions.tmp"

		# Extracting the protein sequences of species mentioned in the species_mapping file. REvolver breaks if the query sequence contains X characters.
		awk -v seqs_dir="$split_seqs_dir" 'BEGIN{RS=">";ORS="";FS=OFS="\n"} FNR==NR{arr[$1]=1;next} $0==""||/#/{next} {species=substr($1,1,5)} arr[species]{protein=$1; $1=""; gsub("\n|X","",$0); print ">"protein,$0,"" > seqs_dir"/"species".fa"; protein="FAIL";species="FAIL"}' "$u_files/sequence_extractions.tmp" "$seqs_file"

		if [[ -f "$u_files/sequence_extractions.tmp" ]]
		then
			rm "$u_files/sequence_extractions.tmp"
		fi
	fi
	
}	# ----------  end of function split_seqs  ----------

split_oma_files ()
{
	### We split the oma-pairs and the oma-seqs file, because they contain too many lines to operate efficiently. Oma-groups does not have so many lines and can be used directly.

	echo "Setting up the used_files directory."

	# Setup the used_files directory.
	pairs_file="$u_files/oma-pairs.txt"
	seqs_file="$u_files/oma-seqs.fa"
	split_pairs_dir="$u_files/pairwise_orth"
	split_seqs_dir="$u_files/splitted_seqs"

	### Splitting the databases into species-categorized parts may be time consuming. Therefore, I made this part interruptable by CNTRL-C.
	# Delete the splitted_pairs directory, if anything fails until the next trap change.
	trap exit_split_pairs SIGINT

	split_pairs

	# Delete the splitted_seqs directory, if anything fails until the next trap change. The previous SIGINT trap function is overridden.
	trap exit_split_seqs SIGINT

	split_seqs

	# Reset the CNTRL-C behavior.
	trap - SIGINT

	# Defuse the exit_splitting function
	unset split_pairs_dir split_seqs_dir

	echo "used_files setup complete."
}	# ----------  end of function split_oma_files  ----------

### Active setup ###
split_oma_files
