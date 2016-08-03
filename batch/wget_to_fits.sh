	#!/bin/bash

proton() {
	directory=http://cta.cppm.in2p3.fr/MeudonPrototype/simtel/simtel_proton_30deg_r1000/
	name=simtel_runmeudon_proton_30tel_30deg_$1
	local_directory=~/Software/outputs/sim_telarray/meudon_proton/
	cd ${local_directory}
#	if [ ! -f ${name} ]; then
#	    wget ${directory}${name} &>${local_directory}${name}_progress.txt
#	fi
	python ~/Software/oxpytools/scripts/write_img_fits -f ${local_directory}${name}.gz -o ${local_directory}${name}.fits &>${local_directory}${name}_progress.txt
}
export -f proton
gamma() {
	directory=http://cta.cppm.in2p3.fr/MeudonPrototype/simtel/simtel_gamma_30deg/
	name=simtel_runmeudon_gamma_30tel_30deg_$1
	local_directory=~/Software/outputs/sim_telarray/meudon_gamma/
	cd ${local_directory}
#	if [ ! -f ${name} ]; then
#	    wget ${directory}${name} &>${local_directory}${name}_progress.txt
#	fi
	python ~/Software/oxpytools/scripts/write_img_fits -f ${local_directory}${name}.gz  -o ${local_directory}${name}.fits  &>${local_directory}${name}_progress.txt
}
export -f gamma

source activate cta
parallel --bar proton ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
parallel --bar gamma ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19