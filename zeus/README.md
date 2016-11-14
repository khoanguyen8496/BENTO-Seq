#	Zeus - Different splicing detector v0.01

Tiny package for detecting different splicing events between 2 conditions

##	Methods and features

-	Index PSI-distribution for all events from all BAM files in current workspace
-	Using Hellinger distance to estimate distance of PSI-distribution among control samples (x) and \
between control and disease samples (y)
-	Assume that no different between control samples and disease samples, x might be equal to y, e.g x~y
-	The distance D of point (x,y) to (d):y=x follow the normal distribution, so the null hypothesis is `H0: d=0`

##	Install

```python
python setup.py install
```

##	Usage
```
zeus.py [-h] [--ThreadsN THREADSN] {index,compute} ...

Zeus - Different exons splicing comparison

positional arguments:
  {index,compute}       Zeus Utilities
    index               Index all splicing results
    compute             Compute different splicing events

optional arguments:
  -h, --help            show this help message and exit
  --ThreadsN THREADSN, -t THREADSN
                        Number of threads
```

###	Indexing

```
usage: zeus.py index [-h] --files FILES [--outfile OUTFILE]

optional arguments:
  -h, --help            show this help message and exit
  --files FILES, -i FILES
                        file contained path of all splicing files in workspace
  --outfile OUTFILE, -o OUTFILE
                        Output binary file of all index splicing result
```
Indexing help accelerate computing different splicing events between conditions in downstream
analysis

-	To index all splicing results from all BAM files:

`zeus.py index --files <input files> -o <output path>`

with `<input files>` could be file paths separated by ","

-	Example

`zeus.py index --files "spliced_a.tsv,spliced_b.tsv,spliced_c.tsv,spliced_d.tsv" -o indexed_db.bin`

###	Computing

```
usage: zeus.py compute [-h] [--outfile OUTFILE] --controls CONTROLS --cases
                       CASES --base-data BASE_DATA --serialized-bin
                       SERIALIZED_BIN

optional arguments:
  -h, --help            show this help message and exit
  --outfile OUTFILE, -o OUTFILE
                        Output binary file of all index splicing result
  --controls CONTROLS, -co CONTROLS
                        Input list of file, separated by "," or file of files
                        controls
  --cases CASES, -ca CASES
                        Input list of file, separated by "," or file of files
                        cases
  --base-data BASE_DATA, -b BASE_DATA
                        Annotation DB as JSON
  --serialized-bin SERIALIZED_BIN, -sb SERIALIZED_BIN
                        serialized bin
```

-	Computing method could be used for detecting different splicing events among 2 conditions

-	Usage
```zeus.py compute --controls <condition 1 files> --cases <condition 2 files> 
        --base-data <path to base_data.json file> --serialized-bin <indexed file> -o <output_tsv>```

with base_data.json is the annotation DB attached with the packages ./zeus/base_data.json

-	Example
```zeus.py compute --controls "spliced_a.tsv,spliced_b.tsv" --cases "spliced_c.tsv,spliced_d.tsv --base_data ./zeus/base_data.json --serialized-bin ./indexed_db.bin -o splicing_result.tsv```

                        
                   
