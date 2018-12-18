# MsPAC
Phase reads, assemble haplotypes and detect SVs

### Installing
```
### Installing MsPAC and it's dependencies
git clone https://github.com/oscarlr/MsPAC.git
cd MsPAC
conda env create -f environment.yml 
conda activate MsPAC
python setup.py install

### Installing cluster package that's needed
cd ..
git clone https://github.com/oscarlr/cluster
cd cluster
python setup.py install
```

### Running the tests
```
cd testing
sh run.sh
```
