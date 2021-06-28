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

### Startup ###

if [ ! -d "./prottrace/" ];then
	echo "Error: Please move to the ProtTrace main directory!"
	exit
fi

work_dir="$(echo $PWD)"
u_files="$work_dir/used_files"

### Minimum version definitions ###

declare -A min_version

min_version['linsi']=6
min_version['hmmscan']=3.1
min_version['hmmpress']=3.1
min_version['blastp']=2.7
min_version['iqtree']=1.6.7.1
min_version['puzzle']=5.3
min_version['Rscript']=3.0
min_version['java']=1.7
min_version['figtree']=1.4

min_version['biopython']=1.78
min_version['DendroPy']=4.5.2

### Function definitions ###

### Anaconda ###

check_conda()
{
	# In an activated conda environment, the string specifying the current
	# environment name is not empty.
	if [ ! -z $CONDA_DEFAULT_ENV ]
	then
		echo "Anaconda environment -- $CONDA_DEFAULT_ENV -- is active."
	else
		echo "No activated Anaconda environment."
	fi
}

install_program ()
{
	if  [ ! -z $CONDA_DEFAULT_ENV ]
	then
		conda install --yes --quiet $1
	else
		python -m pip install $1
	fi
}	# ----------  end of function install_program  ----------


install_R ()
{
	if [ ! -z $CONDA_DEFAULT_ENV ]
	then
		# Install R-base with one-time confirmation. If everything is fine, this one-time confirmation should be enough.
		conda install 'r-base' '-y'
	else
		echo "Please install R."
	fi
}	# ----------  end of function install_R  ----------

### General functions ###

check_program ()
{
	if [ -z "$1"  ] ; then
		echo "$1 could not be found!"
	fi
}	# ----------  end of function check_program  ----------


check_begin ()
{
	echo "Checking $1."
}	# ----------  end of function check_begin  ----------

check_complete ()
{
	echo "Check complete."
}	# ----------  end of function check_complete  ----------

### Elaborate software management functions ###

below_min_version ()
{
	# Accepts strings of numbers separated with dots. Strings 
	# are split by dots to get a series of sub-version numbers.
	# Every sub-version is compared between the min_req and the
	# input v_number, until we find a higher or lower sub-version 
	# number.
	IFS='.' read -ra min_req <<< "$1"
	IFS='.' read -ra v_number <<< "$2"
	
	v_length=${#v_number[@]}
	for (( CNTR=0; CNTR<${#min_req[@]}; CNTR+=1 )); do
	
		# If the minimum requirement ever reaches a higher version depth
		# than the software, the software is below min version.
		if [[ $CNTR -ge $v_length ]] ; then
			return
		fi
		
		# Comparing sub-version numbers.	
		if [[ ${v_number[$CNTR]} -gt ${min_req[$CNTR]} ]] ; then
			false
			return
		fi
		
		if [[ ${v_number[$CNTR]} -lt ${min_req[$CNTR]} ]] ; then
			return
		fi

	done

	# Inconclusiveness means equality. Which is enough for us
	false
	return
}	# ----------  end of function below_min_version  ----------


check_below_version ()
{
	software_name="$1"
	software_version="$2"

	# Warn the user, if the software version is below expectations.
	if below_min_version ${min_version[$software_name]} $software_version; then
		echo "$software_name is below the minimum version of ${min_version[$software_name]}"
		true
	else
		false
	fi

	return
}	# ----------  end of function check_version  ----------

catch_critical_nonident ()
{
	echo "$1 was not found!"
	exit 1
}	# ----------  end of function catch_nonident  ----------

check_names ()
{
	### Checks all passed program names for their existence.

	for prog in "$@"; do
		check_program $prog || echo "$prog was not identified!"
	done
	
}	# ----------  end of function check_names  ----------

### Specific dependency management functions ###

find_python_package ()
{
	software_name=$1

	software_version=$(pip list | grep $software_name | sed 's/  */ /' | cut -d' ' -f2)

	if check_below_version "$software_name" "$software_version"
	then	
		echo "Install necessary!"
		#install_program $software_name
	fi

}	# ----------  end of function find_python_package  ----------

python_dependencies ()
{
	check_begin "for python dependencies."

	prog_name=('biopython' 'DendroPy')

	for p in "${prog_name[@]}"
	do
		find_python_package "$p"
	done

	check_complete

}	# ----------  end of function python_dependencies  ----------

check_terminal_programs ()
{
	### check grep, sed and wget availability
	check_begin "the availability of required terminal functions"
	
	# These programs are checked
	term_prog=('grep' 'sed' 'readlink' 'wget' 'gunzip')

	check_names "${term_prog[@]}"
	
	check_complete

}	# ----------  end of function check_programs  ----------


check_java ()
{
	check_begin "Java"
	# Check whether JAVA is installed for REvolver.
	#check_program $JAVA_HOME

	if type -p java; then
	    echo 'Found java executable in PATH'
	    _java=java
	elif [[ -n "$JAVA_HOME" ]] && [[ -x "$JAVA_HOME/bin/java" ]];  then
	    echo 'Found java executable in JAVA_HOME'     
	    _java="$JAVA_HOME/bin/java"
	else
	    echo 'Java was not found!'
	fi
	
	if [[ "$_java" ]]; then
	    version=$("$_java" -version 2>&1 | awk -F '"' '/version/{print $2}')
	    if below_min_version "${min_version["java"]}" "$version" ; then
	        echo "Java version is lower than java 1.7 / openjdk 7.0!"
	    fi
	fi

	check_complete

}	# ----------  end of function check_java  ----------


resolve_alias ()
{
	### Legacy function. Currently unused. 
	grep "$1" ~/.bash_aliases | cut -d'=' -f2 

}	# ----------  end of function resolve_alias  ----------

find_software ()
{
	### Tries to check the existence and the version of the requested software. If it does not exist yet, the function tries to install it.

	software_name=$1

	check_begin $software_name
	
	### If the software is set as an alias, we have to dig the path. ###
	alias_path=$(resolve_alias $software_name)
	
	### If the software name represents a file, we can extract the path from type.
	if [[ $(type -t $software_name) == "file" ]]
	then
		file_path=$(type -p $software_name)
	fi
	software_path="${alias_path:-$file_path}"

	### Retrieve the software version. ###
	software_version=""

	# hmmscan and treepuzzle, for example, do not provide a --version output. Instead, we must read the first match on the -h induced help page.
	tried_arguments=('--version' '-h')
	for arg in "${tried_arguments[@]}"
	do
		# Some programs may show version information with stderr. This is why we combine stderr with stdin with 2>&1.
		software_version=$($software_path $arg 2>&1 | grep -E -m1 -o "[0-9]+\.[0-9]+\.?[0-9]*\.?[0-9]*" | head -n1)
		if [ ! -z $software_version ]
		then
			break
		fi
	done


	if check_below_version "$software_name" "$software_version"
	then
		if [ "$software_name" == "Rscript" ]
		then
			echo "Please install r-base!"
			# Anaconda can take long with resolving the environment. This piece of code is kept in at leasure. It should work by just uncommenting the line below.
			#install_R
		fi
	fi
	
	check_complete "$software_name"
	echo "-----------------------------------------------------------------------"
	
	# We have to use # as sed command delimiters, because software_path contains backslashes.
	sed -ri 's#^('"${software_name}"':).*#\1'"${software_path}"'#' prog.config
	
	# We need to unset the variables, or otherwise they bleed to those we cannot find.
	unset software_name software_path software_version alias_path file_path
}	# ----------  end of function find_software  ----------

check_config_programs ()
{
	### Checks immediate access to programs in order to add them to the config file.
	### Programs missing here can be added manually later.

	check_begin "for remaining third-party software. May correct above shown errors about missing programs"

	prog_name=('hmmscan' 'hmmpress' 'blastp' 'iqtree' 'puzzle' 'linsi' 'figtree' 'Rscript')

	for p in "${prog_name[@]}"
	do
		find_software "$p"
	done

	check_complete	

}	# ----------  end of function check_config_programs  ----------


download_oma_db ()
{
	pairs_link='https://omabrowser.org/All/oma-pairs.txt.gz'
	groups_link='https://omabrowser.org/All/oma-groups.txt.gz'
	seq_link='https://omabrowser.org/All/oma-seqs.fa.gz'

	wget --no-verbose -N -O "$u_files/oma-pairs.txt.gz" "$pairs_link"
	wget --no-verbose -N -O "$u_files/oma-groups.txt.gz" "$groups_link"
	wget --no-verbose -N -O "$u_files/oma-seqs.fa.gz" "$seq_link"
	gunzip "$u_files/oma-pairs.txt.gz" "$u_files/oma-groups.txt.gz" "$u_files/oma-seqs.fa.gz"

}	# ----------  end of function download_oma_db  ----------


download_pfam_db ()
{
	pfam_link='ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'

	wget -N -O "$u_files/Pfam-A.hmm.gz" "$pfam_link"
	gunzip "$u_files/Pfam-A.hmm.gz"
	hmmpress "$u_files/Pfam-A.hmm"

}	# ----------  end of function download_pfam_db  ----------

create_config ()
{
	if [ ! -f "$u_files/oma-groups.txt" ] && [ ! -f "$u_files/oma-seqs.fa" ]
	then
		download_oma_db
	fi

	if [ ! -f "$u_files/Pfam-A.hmm" ]
	then
		download_pfam_db
	fi

	perl ./prottrace/setup/create_conf.pl -quiet -name "prog.config"

}	# ----------  end of function create_config  ----------


split_oma ()
{
	# Split the OMA database into smaller pieces for much faster lookup of individual species.The script has been outsourced to make separate reruns possible.
	"$work_dir/prottrace/setup/split_oma.sh"

	# Change the oma database paths to the new directories, if they exist.
	if [[ -d "$u_files/splitted_seqs" ]]
	then
		split_seqs_path="$u_files/splitted_seqs"
		sed -ri 's#(path_oma_seqs:).*#\1'"$split_seqs_path"'#' prog.config
	fi

	if [[ -d "$u_files/pairwise_orth" ]]
	then
		split_pairs_path="$u_files/pairwise_orth"
		sed -ri 's#(path_oma_pairs:).*#\1'"$split_pairs_path"'#' prog.config
	fi


}	# ----------  end of function split_oma  ----------

### Executing tests ###

check_conda
python_dependencies
check_terminal_programs
check_java

### Creating default directories ###

[ -d "output" ] || mkdir "output"
[ -d "cache" ] || mkdir "cache"
[ -d "distances" ] || mkdir "distances"

### Executing the older script that generates the config file. It remains primarily useful to download databases. ###

create_config

# Fills in missing programs
check_config_programs

### Active setup ###
split_oma
