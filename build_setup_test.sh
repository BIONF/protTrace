#!/bin/bash

#===============================================================================
#
#          FILE: build_test.sh
#
#         USAGE: ./build_test.sh
#
#   DESCRIPTION: Automatically builds, reinstalls and runs ProtTrace.
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Dominik Perisa (DP), dominik.perisa@bio.gmail.com
#  ORGANIZATION: Goethe University Frankfurt
#       CREATED: 04/15/2021 07:42:53 PM
#      REVISION:  ---
#===============================================================================

# Builds the package. A wheel is added to the ./dist directory.
python -m build

# Retrieves the latest wheel.
unset -v latest
for file in $(find dist/ -maxdepth 1 -name ProtTrace-*-py3-none-any.whl); do
	[[ $file -nt $latest ]] && latest=$file
done

# Installs the wheel into the current (conda) environment.
python -m pip install --force-reinstall "${latest}"

# Runs the program.
prottrace -i "example/test.id" -c "prog.config"
# Fdog is rather tested with the fasta option.
#prottrace -f "example/test.fasta" -c "prog.config"
