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

work_dir="$(echo $PWD)"
setup_pairwise=false

# Get arguments.
while getopts w:p opt; do
	case ${opt} in
		p ) setup_pairwise=true;;
	esac
done

check_program ()
{
	if [ -z "$1"  ] ; then
		echo "$1 could not be found!"
	fi
}	# ----------  end of function check_program  ----------


install_program ()
{
	if  [ $1 == "conda" ]; then
		conda install -c conda-forge $2
	else
		echo "The first argument must specify the install method!"
	fi
}	# ----------  end of function install_program  ----------

check_begin ()
{
	echo "Checking $1."
}	# ----------  end of function check_begin  ----------

check_complete ()
{
	echo "Check complete."
}	# ----------  end of function check_complete  ----------


below_min_version ()
{
	# Accepts strings of numbers separated with dots. Strings 
	# are split by dots to get a series of sub-version numbers.
	# Every sub-version is compared between the min_req and the
	# input v_number, until we find a higher or lower sub-version 
	# number.
	IFS='.' read -ra min_req <<< $1
	IFS='.' read -ra v_number <<< $2

	v_length=${#v_number[@]}
	for (( CNTR=0; CNTR<${#min_req[@]}; CNTR+=1 )); do
	
		# Assuming the minimum requirement to be less or equally 
		# specific, a dragged on comparison means equal sub-version
		# numbers. Which is sufficient to falsify the need to warn
		# the user.
		if [[ $CNTR -ge $v_length ]] ; then
			false
		fi
		
		# Comparing sub-version numbers.	
		if [[ ${v_number[$CNTR]} -gt ${min_req[$CNTR]} ]] ; then
			false
		fi
		
		if [[ ${v_number[$CNTR]} -lt ${min_req[$CNTR]} ]] ; then
			return
		fi

		# Inconclusiveness means equality.
	done
}	# ----------  end of function min_version  ----------

check_programs ()
{
	
	catch_nonident ()
	{
		echo "$1 was not found!"
		exit 1
	}	# ----------  end of function catch_nonident  ----------

	### check grep, sed and wget availability
	check_begin "the availability of required terminal functions"
	
	# These programs are checked
	nextprog=('grep' 'sed' 'readlink' 'wget')

	for prog in nextprog; do
		check_program $prog || echo "$prog was not identified!"
	done
	
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
	    if ! min_version "7.0" "$version" ; then
	        echo "Java version is lower than java 1.7 / openjdk 7.0!"
	    fi
	fi

	check_complete

}	# ----------  end of function check_java  ----------

check_programs
check_java

echo "Setting up the used_files directory."
# Setup the used_files directory.
u_files="$work_dir/used_files"

if [ "$setup_pairwise" == true ] ; then
	if [ -f $u_files/oma_pairs.txt ];then
		pairwise_spec_dir="$u_files/pariwise_orth"

		if [ -d $pairwise_spec_dir ];then
			mkdir $pairwise_spec_dir
		fi

		awk -v pairwise_dir=$pairwise_spec_dir 'BEGIN{FS=OFS="\t"} {s1=substr($1,1,5); s2=substr($2,1,5); print $1,$2,s1,s2 > pairwise_dir"/"s1"_"s2".txt"}' "$u_files/oma_pairs.txt"
	fi
fi
echo "used_files setup complete."
