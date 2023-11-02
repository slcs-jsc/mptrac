#! /bin/bash

# Create directories...
rm -rf data && mkdir -p data

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH

# Path to log files...
log=./data

# Pathway to the application...
atm_conv=../../src/atm_conv

# Path to the atmosphere data that containts air parcel positions...
atm_in=./data.ref/atm_input.tab

# Path to the output atmosphere data...
atm_out_nc=./data/atm_output.nc

# Path to control file...
ctl=./data/trac.ctl

# Create control parameter file...
cat > ./data/trac.ctl <<EOF
EOF

# Run the application...
echo "Test reading atm_type 0 data and write it to atm_type 4 data."
echo ${atm_conv} ${ctl} ${atm_in} 0 ${atm_out_nc} 4
${atm_conv} ${ctl} ${atm_in} 0 ${atm_out_nc} 4 >> ${log}/log_0to4.txt

# Path to the atmosphere data that containts air parcel positions...
atm_in=${atm_out_nc}

# Path to the output atmosphere data...
atm_out_tab=./data/atm_output.tab

# Path to the reference...
atm_ref_tab=./data.ref/atm_output.tab

# Path to control file...
ctl=./data/trac.ctl

# Create control parameter file...
cat > ./data/trac.ctl <<EOF
EOF

# Run the application...
echo "Reading atmospheric data of atm_type 4 and write it as atm_type 0."
echo ${atm_conv} ${ctl} ${atm_in} 4 ${atm_out_tab} 0
${atm_conv} ${ctl} ${atm_in} 4 ${atm_out_tab} 0 > ${log}/log_4to0.txt

# Compare files...
error=0
diff -q -s ${atm_out_tab} ${atm_ref_tab}

if [ $? -ne 0 ]; then
	error=1
	cat ${log}/log_0to4.txt
	cat ${log}/log_4to0.txt
	exit $error
fi

#rm ${atm_out_tab}
rm ${atm_out_nc}

rm ${ctl}
rm ${log}/*.txt

exit $error
