The inputs folder contains a sample parsed reference template for chromosome 21 and a sample parsed reads dataset simulated with SimLord (150bp long, 100K reads; error rate 0.01) that can be used for testing the prototype. 

Note: The preprocessing run will generate several auxiliary files for reference index. Additionally, User is expected to specify the path of the main directory in preprocessing.py, alignment.py, and postprocessing.py files.

options for preprocessing, alignment, and postprocessing:

```
--chrnum: chromosome number 
--nr: total number of reads 
--rl: read length 
--nrg: number of read groups 
--nrb: number of read batches 
--tcn: template chunk number 
--rbn: read batch number 
```

Required python modules:

```
pip install numpy 
pip install mpi4py
```
```
cd codes 
```

Example commands:

```
python3 preprocessing.py --chrnum 21 --nr 100000 --rl 150 --nrg 2 --nrb 8
python3 alignment.py --chrnum 21 --rl 150 --nrg 2 --nrb 8 --tcn 0 --rbn 0
python3 postprocessing.py --chrnum 21 --rl 150 --nrg 2 --nrb 8 --tcn 0 --rbn 0
```


This code is a prototype and is not intended for large-scale commercial applications.
