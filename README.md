# Motif_Finding_DNA
CS 412 Final Project. UIUC. Spring 2020
Team:
Gowtham Kuntumalla
Devin Tark
Xaiowei Lyu

Instructions:
a. The runtime for this program is long - O(hours). 
b. Instructions are for MacOS zsh terminal. Linux Bash commands are similar.

1. make the py files excutable
chmod 755 *py 

2. check status:
stat -f "%OLp" benchmark.py
stat -f "%OLp" motif_finder.py
stat -f "%OLp" evaluate.py

3. run the benchmarks program using
nohup python3 benchmark.py > benchmark_output.log &

3. run the EM algo & evaluation program using
nohup python3 evaluate.py > evaluate_output.log &

4. check status using 
ps | grep [PID]

if you want to exit the process:
kill [PID]
