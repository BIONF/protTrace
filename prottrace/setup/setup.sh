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

### check grep, sed and wget availability
check_begin "the availability of required terminal functions"
grepprog='grep'
check_program $grepprog
sedprog='sed'
check_program $sedprog
readlinkprog='readlink'
check_program $readlinkprog
wgetprog='wget'
check_program $wgetprog

check_complete

check_begin "non-python/bash requirements"
# Check whether JAVA is installed for REvolver.
check_program $JAVA_HOME

check_complete

echo "Setting up the used_files directory."
# Setup the used_files directory.
u_files="$work_dir/used_files"

if [ "$setup_pairwise" == true ] ; then
	if [ -f $u_files/oma_pairs.txt ];then
		pairwise_spec_dir="$u_files/pariwise_orth"

		if [ -d $pairwise_spec_dir ];then
			mkdir $pairwise_spec_dir
		fi

		#awk -v pairwise_dir=$pairwise_spec_dir 'BEGIN{FS=OFS="\t"} {s1=substr($1,1,5); s2=substr($2,1,5); print $1,$2,s1,s2 > pairwise_dir"/"s1"_"s2".txt"}' "$u_files/oma_pairs.txt"
	fi
fi
echo "used_files setup complete."
