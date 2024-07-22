# Sameer
Program for getting mapping statistics, generate consensus and to get low frequency variats from a sam file.

## Installation
```
git clone https://github.com/vbsreenu/Sameer.git
cd Sameer
cc sameer_new.c -o SAMEER_NEW
chmod +x SameerReport
chmod +x *.R 
```

** Set the path to the installation directory to access these programs. (This is dependent on your shell) **


## Usage
```
SameerReport input.sam annotation.tsv Ref.fa
```


## Demo run

```
cd demo
unzip 1-MRU25010-M-mapped.sam.zip
cd ..
./SameerReport demo/1-MRU25010-M-mapped.sam demo/MRU25010-M.anno demo/MRU25010-M.fa
```
Here 
1st input is a SAM file
2nd input is an annotation file. This is a tab seperated file with four columns, ORF start, ORF end, ORF orientation and Name.
3rd input is the reference file used for mapping


## Output
This will produce a PDF file with the detailed summary and text file of amino acid variations (`variations.txt`). 
In the demo output it should look like...

```
Glycoprotein    21-3614 1198
1 M M   M: 99.66(583)
2 Y Y   Y: 99.37(634)
3 V V   V: 100.00(1120)
4 L L   L: 99.93(1399)
5 L L   L: 99.93(1422)
6 T T   T: 100.00(1498)
7 I I   I: 99.93(1530)
8 L L   L: 99.94(1687)
9 I I   I: 99.88(1724)
10 S S  S: 99.89(1895)
```

PDF output can be see in demo direcory.





## Dependencies:
- Latex
- R
- gcc (or cc)

