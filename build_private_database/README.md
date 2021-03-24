### usage
1. generate the private database
```
./run_genSubsetDB.sh
```

2. obtain target haplotypes in the private database
```
python genSubsetDB_getHaplotype.py
```

3. obtain haplotypes not in the private database to build the null distribution
```
python genSubsetDB_getHaplotype_null.py
```
