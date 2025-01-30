# The S1A serine proteases 
echo "S1A serine protease Calculations:" > Outputs/s1A_halabi.log
./scaProcessMSA.py Inputs/s1Ahalabi_1470_nosnakes.an -s 3TGI -c E -t -n >> Outputs/s1A_halabi.log 
./scaCore.py Outputs/s1Ahalabi_1470_nosnakes.db >> Outputs/s1A_halabi.log 
./scaSectorID.py Outputs/s1Ahalabi_1470_nosnakes.db >> Outputs/s1A_halabi.log 
# Beta-lactamase 
echo "Beta-lactamase Calculations:" > Outputs/PF13354.log
./scaProcessMSA.py Inputs/PF13354_full.an -s 1FQG -c A -f 'Escherichia coli' -t -n >> Outputs/PF13354.log
./scaCore.py Outputs/PF13354_full.db >> Outputs/PF13354.log
./scaSectorID.py Outputs/PF13354_full.db >> Outputs/PF13354.log
# G-protein - this analysis is run with two alignments - the full PFAM 
# alignment (PF00071_full) and the PFAM alignment filtered to remove several 
# N-terminal truncation mutants. PF00071_rd2 is the aligment discussed in the
# manuscript.
echo "G-protein calculations:" > Outputs/PF00071.log
./scaProcessMSA.py Inputs/PF00071_full.an -s 5P21 -c A -f 'Homo sapiens' -t -n >> Outputs/PF00071.log
./scaCore.py Outputs/PF00071_full.db >> Outputs/PF00071.log
./scaSectorID.py Outputs/PF00071_full.db >> Outputs/PF00071.log
echo "G-protein calculations:" > Outputs/PF00071_rd2.log
./scaProcessMSA.py Inputs/PF00071_rd2.an -s 5P21 -c A -f 'Homo sapiens' -t -n >> Outputs/PF00071_rd2.log
./scaCore.py Outputs/PF00071_rd2.db >> Outputs/PF00071_rd2.log
./scaSectorID.py Outputs/PF00071_rd2.db >> Outputs/PF00071_rd2.log
# DHFR - this analysis is also run with two alignments for comparison - 
# the full PFAM alignment (PF00186_full.an) and a manually curated alignment 
# (DHFR_PEPM3.an)  
echo "DHFR Calculations:" > Outputs/PF00186.log
./scaProcessMSA.py Inputs/PF00186_full.an -s 1RX2 -c A -f 'Escherichia coli' -t -n >> Outputs/PF00186.log
./scaCore.py Outputs/PF00186_full.db >> Outputs/PF00186.log
./scaSectorID.py Outputs/PF00186_full.db >> Outputs/PF00186.log
echo "DHFR Calculations:" > Outputs/DHFR_PEPM3.log
./scaProcessMSA.py Inputs/DHFR_PEPM3.an -s 1RX2 -c A -t -n >> Outputs/DHFR_PEPM3.log
./scaCore.py Outputs/DHFR_PEPM3.db >> Outputs/DHFR_PEPM3.log
./scaSectorID.py Outputs/DHFR_PEPM3.db >> Outputs/DHFR_PEPM3.log

