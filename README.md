The inputs folder contains a sample parsed reference template for chromosome 21 and a sample parsed reads dataset simulated with SimLord (150bp long, 100K reads; error rate 0.01) that can be used for testing the prototype. Note that the preprocessing run will generate several auxiliary files for reference index. 

Execute the following commands to create and activate virtual environment and install the required packages:
```
python3 -m venv venv  
source venv/bin/activate
pip install -r requirements.txt
```

Options for preprocessing, alignment, and postprocessing:

```
--chrnum: chromosome number 
--nr: total number of reads 
--rl: read length 
--nrg: number of read groups 
--nrb: number of read batches 
--tcn: template chunk number 
--rbn: read batch number 
```


```
cd codes 
```

Example commands for testing your installation:

```
python3 preprocessing.py --chrnum 21 --nr 100000 --rl 150 --nrg 2 --nrb 100
mpirun -np 2 python3 alignment.py --chrnum 21 --rl 150 --nrg 2 --nrb 100 --tcn 0 --rbn 0
python3 postprocessing.py --chrnum 21 --nr 100000 --rl 150 --nrg 2 --nrb 100 --tcn 0 --rbn 0
```
The job size for the above alignment task is a template size of 1M and 1000 reads. It takes about 2.5 minutes to run on an 8-cores 8 GB RAM MacBook M1 machine.

This code is a prototype and is not intended for large-scale commercial applications.
