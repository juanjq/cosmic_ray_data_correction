#! /bin/bash

file="/fefs/aswg/workspace/juan.jimenez/lst1_systematics/bash_energy_ranges/jobs_information.txt"

# read the lines of the given file
while read -r line; do

# only operating if the line is not commented with #
if  [[ "${line:0:1}" != '#' ]]; then
energy_str="$line"

echo "#! /bin/bash" > compute_lc_tmpjob.sh
echo "python /fefs/aswg/workspace/juan.jimenez/lst1_systematics/bash_energy_ranges/script_compute_lc.py compute '$energy_str'" >> compute_lc_tmpjob.sh

echo -e "Sending energy range $energy_str to compute the LC"
sbatch -p short --output="/fefs/aswg/workspace/juan.jimenez/lst1_systematics/bash_energy_ranges/output_slurm/slurm-%j.out" compute_lc_tmpjob.sh

fi

done < $file
