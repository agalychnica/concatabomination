#/bin/bash
export DESTINATION_PATH=~/bin/;

if [ ! -e $DESTINATION_PATH ]; then 
	echo "ERROR: $DESTINATION_PATH does not exist, please create this folder or change the definition on line 2 of this installation script"; 
	exit 1; 
fi;

# COMPILE and install concatabomb
cc -o concatabomb concatabomb_v1.2.c
if [ ! -e concatabomb ]; then 
	echo "ERROR: could not sucessfully compile concatabomb; please check gcc is installed and on your PATH"; 
	exit 1;
else 
	cp concatabomb $DESTINATION_PATH;
fi 

# COMPILE and install pairwise_incompatibility
cc -o pairwise_incompatibility pairwise_incompatibility.c
if [ ! -e pairwise_incompatibility ]; then 
	echo "ERROR: could not sucessfully compile pairwise_incompatibility; please check gcc is installed and on your PATH"; 
	exit 1;
else 
	cp pairwise_incompatibility $DESTINATION_PATH;
fi 

# Install PerlEQ
if [ ! -e PerlEQ_orig_version.pl ]; then
	echo "ERROR: PerlEQ_orig_version.pl is missing from the installation, please download a new version of the pipline";
	exit 1;
else
	cp PerlEQ_orig_version.pl $DESTINATION_PATH;
fi

# Install COMPASS
if [ ! -e COMPASS.tar ]; then
	echo "ERROR: COMPASS is missing from the installation, please download a new version of the pipline";
	exit 1;
else
	tar -xvf COMPASS.tar;
	cd COMPASS/src
	rm COMPASS
	make;
	if [ ! -e COMPASS ]; then
		echo "ERROR: There was a problem installing COMPASS; please see the installation guide in the folder \"COMPASS\" to see if your system meets its requirements"
		exit 1; 
	else
		cp COMPASS $DESTINATION_PATH;
		cd ../../;
	fi
fi

# INSTALL concatabomination_pipeline_v4_parallel

if [ ! -e concatabomination_pipeline_v4.1_parallel ]; then
	echo "ERROR: concatabomination_pipeline_v4_parallel is missing from the installation, please download a new version of the pipline";
	exit 1;
else
	chmod a+x concatabomination_pipeline_v4.1_parallel;
	cp concatabomination_pipeline_v4.1_parallel $DESTINATION_PATH;
fi


cat << EOF


Success! 
Test the installation of the pipline by typing the following command:

	concatabomination_pipeline_v4.1_parallel -i Gauthier_1986.coding.nex

Type: 
	concatabomination_pipeline_v4.1_parallel -h
to see all the possible options.

EOF
exit 0;





