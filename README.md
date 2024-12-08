# Colormap-Project
The goal of this project is to re implement the methodology presented in the paper "CoLoRMap: Correcting Long Reads by Mapping short reads" by Ehsan Haghshenas, Faraz Hach, S. Cenk Sahinalp and Cedric Chauve


# Dependancies

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

working on this section