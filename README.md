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
```sudo apt-get install libboost-graph-dev```

# Usage

## Snakemake pipeline

A folder with 4 files
1. ```<short reads 1>.fastq```
2. ```<short_reads_1>.fastq```
3. ```<long reads>.fasta```
3. ```<reference genome>.fasta```

parameterizes the Snakefile. These values can be specified in the snakefile.


## ```colormap.cpp```
This file can be used directly to correct long reads. It takes 2 command line arguments:
1. ```<long_reads>.fasta```  
this is just the relative path to the long reads which are 
2. a "raw alignment file"



