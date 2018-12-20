# MsPAC
Phase reads, assemble haplotypes and detect SVs

### Usage
```
MsPAC phase-bam run.cfg
MsPAC prep-reads run.cfg
MsPAC assembly run.cfg
MsPAC sv-calling run.cfg
```

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
git clone https://github.com/oscarlr/cluster.git
cd cluster
python setup.py install
```

### Running the tests
```
cd testing
sh run.sh
```
