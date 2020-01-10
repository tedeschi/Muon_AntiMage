#!/bin/bash

#remove the old stuff
rm Run-*sh
rm Run-*mac


for run in {1..10}
do
	echo $run
	macroname="Run-$run.mac"
	outputfilename="Run-$run"
	shellname="Run-$run.sh"

  #write the macro geant4 will call
  echo "/file/setFileName $outputfilename" > $macroname
	echo "/run/beamOn 10000" >> $macroname

  #write the script to be submitted via slurm
  echo "#!/bin/bash -l" > $shellname
 	echo " " >> $shellname
 	echo "#SBATCH -p debug" >> $shellname
 	echo "#SBATCH -t 48:00:00" >> $shellname
 	echo "#SBATCH -J 2vbbAlpha-$detector-$run" >> $shellname
  echo "#SBATCH -N 1" >> $shellname
 	echo "./Muon_AntiMage $macroname" >> $shellname
 	
 	#submit
	sbatch -p debug -t 48:00:00 $shellname &

  #give it a few seconds before submitting next script so random numbers are different
	sleep 2
done

