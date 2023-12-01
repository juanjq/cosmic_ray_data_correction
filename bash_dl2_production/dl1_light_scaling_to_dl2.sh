#! /bin/bash

file="/fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/bash_dl2_production/jobconfig.txt"

# read the lines of the given file
while read -r line; do

# only operating if the line is not commented with #
if  [[ "${line:0:1}" != '#' ]]; then
runsubrun_str="$line"

echo "#! /bin/bash" > dl1_light_scaling_to_dl2_tmpjob.sh
echo "python /fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/bash_dl2_production/script_process_dl1a_to_dl2_run.py '$runsubrun_str'" >> dl1_light_scaling_to_dl2_tmpjob.sh


# Split the string into an array using - as the delimiter
# splitted='-' read -ra subrun_str <<< "$subrun_str"
IFS='-' read -r -a splitted <<< "$runsubrun_str"
# Get the first element from the array
n_hours="${splitted[0]}"

echo -e "Sending job $runsubrun_str to convert dl1 to dl2 with different light scaling"
sbatch -p short --mem=15000 --output="/fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/bash_dl2_production/output_slurm/slurm-%j.out" --begin="now+${n_hours}hour" dl1_light_scaling_to_dl2_tmpjob.sh

fi

done < $file
