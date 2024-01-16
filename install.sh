#!/bin/bash

### THE CURRENT PATH
mstoolpath=$(pwd)


### MESSAGES TO USERS
printf '\n\n\n...... Thank you for using mstool\n'
echo   '...... Adding mstool to the environment (~/.bashrc and ~/.zshrc)'
echo   "export mstool=$mstoolpath/"
echo   "export PYTHONPATH=$mstoolpath/src:\$PYTHONPATH"
printf '\n\n\n'


########## CMDLINE ##########
export mstool=$mstoolpath/
export PYTHONPATH=$mstoolpath/src:$PYTHONPATH

echo   "export mstool=$mstoolpath/" >> ~/.bashrc
echo   "export PYTHONPATH=$mstoolpath/src:\$PYTHONPATH" >> ~/.bashrc
echo   "export mstool=$mstoolpath/" >> ~/.zshrc
echo   "export PYTHONPATH=$mstoolpath/src:\$PYTHONPATH" >> ~/.zshrc
#############################


### MESSAGES TO USERS
echo   "...... mstool requires the following python packages:"
printf "openmm, sqlite3, cython, numpy, pandas\n\n\n\n"
echo   "...... Testing whether you have these packages"


### TESTING PYTHON PACKAGES
for package in openmm sqlite3 cython numpy pandas; do
    if python -c "import $package" &> /dev/null; then
        echo "$package is found"
    else
        echo "$package is NOT found"
    fi
done
printf "\n\n\n"


### MESSAGES TO USERS
echo "...... mstool uses cython to calculate a distance matrix"
echo "...... Building a cython distance matrix program"


########## CMDLINE ##########
cd src/mstool/lib
python setup.py build_ext --inplace &> /dev/null
#############################


### MESSAGES TO USERS
#if test -e *.so || test -e *.pyd; then
count=`ls -1 *.so 2>/dev/null | wc -l`
if [ $count != 0 ]; then
    echo "Identifying src/mstool/lib/*.so"
    echo "distance.distance_matrix is successfully installed"
else
    echo "distance.distance_matrix is NOT successfully installed"
    echo "You need to install this to make a martini system"
    echo "Otherwise, it is fine without this"
fi 


########## CMDLINE ##########
cd ../../../
#############################


printf '\n\n\n...... All done! Thank you for using mstool\n\n\n\n'

