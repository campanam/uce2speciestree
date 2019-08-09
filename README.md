# uce2speciestree
Pipeline script to generate input files for ASTRAL and SVDquartets from phyluce nexus alignments  

Michael G. Campana, 2019  
Smithsonian Conservation Biology Institute  
Contact: <campanam@si.edu>  

Pipeline script to generate input files for [ASTRAL](https://github.com/smirarab/ASTRAL) and [SVDquartets + PAUP*](https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/) from a collection of [phyluce](https://phyluce.readthedocs.io/en/latest/) UCE loci alignments in nexus format.  

## License  
The software is made available under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Installation  
In the terminal:  
`git clone https://github.com/campanam/uce2speciestree`  
`cd uce2speciestree; chmod +x uce2speciestree.rb`  

## Configuration  
This pipeline is hard-coded for the Smithsonian High Performance Computing Cluster ('Hydra') using SGE. It requires RAxML and Ruby and loads those software as modules. The code will not work on other clusters without manual editing of the script. Assuming a similarly configured SGE cluster, the modules and executables to be changed are located in the write_array_qsub method. The software is provided 'as-is' as a useful reference, but there are no current plans to develop it for general use.    

## Execution  
`ruby uce2speciestree.rb` will show the help menu for execution.  

## Citation  
Campana, M.G. 2019. uce2speciestree. https://github.com/campanam/uce2speciestree.
