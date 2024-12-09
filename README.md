# Colormap-Project
The goal of this project is to re implement the methodology presented in the paper "CoLoRMap: Correcting Long Reads by Mapping short reads" by Ehsan Haghshenas, Faraz Hach, S. Cenk Sahinalp and Cedric Chauve


# Dependencies

- Entrez Direct: used for downloading reference sequence NC 000913  
    ```apt install ncbi-entrez-direct```

- BWA - see here for installation : https://github.com/lh3/bwa. This code base expects it to be in ```$PATH```  

- zlib.h  
```apt-get install zlib1g-dev```

- samtools  
```apt install samtools```

- Boost: C++ library used in the codebase to for store graphs,find connected components,run dijstra's shortest path etc.  
```apt-get install libboost-graph-dev```

- BLASR - used to align long reads to reference genome  
```apt install -y blasr```

# Usage

## Snakemake pipeline

#### Input:
The main 4 parameter for ```Snakefile``` are 

1. ```folder```: The name of the target folder which contains (ex. test_data/)
2. ```<short reads 1>.fastq```: The name of one of fastq files in ```folder``` (ex. ```ill_1.fastq```)
3. ```<short_reads_1>.fastq```: The name of the other fastq file in ```folder``` (ex. ```ill_2.fastq```)
4. ```<long reads>.fasta```: The name of the fasta file in ```folder``` (ex. ```pac.fasta```)

The other parameters are
-  ```test_name```: The suffix of the file containing the corrected long reads. More spefically, the corrected long reads will be stored in ```<folder>/lr_corr_<test_name>.fasta```. This is an id that is intended to be used to distiguish the output files produced as the ```colormap.cpp``` is adjusted.
-  ```correct_singletons```: When set to ```"no"```, then a short read $s$ which has been mapped to a long read $l$ and is not adjacent to any other short reads mapped to $l$ will not be used to correct $l$. Other wise such short reads will be used to correct $l$

#### Output:

This pipeline produces the file 


## ```colormap.cpp```
This file can be used directly to correct long reads. It takes 2 command line arguments:
1. ```<long_reads>.fasta```  
this is just the relative path to the long reads which are 
2. a "raw alignment file"



