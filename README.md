# ORFtranslator
---
ORF prediction from DNA or RNA sequence  

### Requirement
Python 3
Install Biopython


### Installation requirement on Conda environment
```
conda install -c conda-forge biopython
```
[Biopython Install guide with Anaconda](https://anaconda.org/conda-forge/biopython)
  
### Usage with default parameters
```
python ORFtranslator.py -i <input fasta file> -o <output file name> 
```

### Usage with parameters
```
python ORFtranslator.py -i <input fasta file> -o <output file name> -min <min length of AA> --start < Allow start codons, 0:ATG, 1:near-cognate codons, 2: Any codons> --codontable <Genetic code (1-31)> --strand < plus|minus|both >
```

### Example usages
min length 10aa  
start AUG  
+ strand only
Standard Codon Tables
```
python ORFtranslator.py -i in.fa -o out -min 10 --start 0 --codontable 1 --strand plus
```

min length 30aa  
start with near-cognate codons
reverse strand only
Standard Codon Tables
```
python ORFtranslator.py -i in.fa -o out -min 30 --start 1 --codontable 1 --strand minus
```

If you want selected near-cognate codons (non-AUG start codons), fix the file "Initiation_Codon_lists.txt"  

Codon table: see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for detail




