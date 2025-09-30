#!/bin/bash

python test_meep_socketmolecule_tlsrelaxation.py &

echo "After launching MEEP, sleeping 20 seconds to wait for the socket server to start..."
sleep 20

mxl_driver.py --model tls --port 31886 --param "omega=0.242, mu12=187, orientation=2, pe_initial=0.01" --verbose