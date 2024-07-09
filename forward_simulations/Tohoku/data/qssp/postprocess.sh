#!/bin/bash 


# NOTE: The python script qssp_to_sac_loc requires obspy 
# For me the easiest way to do this is to use a conda environment that 
# I know has it in. YOu will need to set up your own conda.sh path 
# Use command conda info | grep -i 'base environment'
# to get the path and then get the etc/profile.... 
# If not you can just run the python script however you normally would
# separately, then run the rest of this bash script after. 
source /usr/local/Caskroom/miniconda/base/etc/profile.d/conda.sh
conda activate SPECFEMX_BENCHMARKS

master_path=$(pwd)
echo $master_path
dir_path="$master_path/prem_elastic/"
ZACC="$dir_path/Zacc"
GRAV="$dir_path/Grav"

# Sac executable: 
SACHOME='/Users/eaton/sac'
sac="$SACHOME/bin/sac"

# Copy the qssp_to_sac python script 
cp ../../../../qssp_to_sac.py qssp_to_sac_loc.py

# Correct the inputs for this simulation: 
sed -i '' "s:^qssp_dir.*:qssp_dir = '$dir_path':" qssp_to_sac_loc.py


# Check if directories exist and if they do then clean them, if not make them 
if [ -d $ZACC ]; then
  echo "Zacc directory exists."
  #rm $ZACC/* 
else 
  mkdir $ZACC
fi
if [ -d $GRAV ]; then
  echo "Grav directory exists."
  #rm $GRAV/* 
else 
  mkdir $GRAV
fi

# Now convert to sac using python code: 
python3 -m qssp_to_sac_loc

# Apply the SAC macro to process the sac files
cd $GRAV 
$sac ../macro_qssp_sac
rm *_tmp.SAC

cd $ZACC
$sac ../macro_qssp_sac
rm *_tmp.SAC