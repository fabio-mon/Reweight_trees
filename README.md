# HH reweight tool     
Tool to reweight any input benchmark to any output benchmark.   

## reweight_HH.py     
Contains the HHreweighter class which is the core of the reweight tool. 
In order to construct a HHreweighter object, the following data are required:
* the input node 
* the requested output node      
* the path of the input mHH shapes
` myreweighter = reweight_HH.HHreweighter(inputnode, outputnode, inputshapesNLO, inputshapesLO, inputshapesLOfake) `      
the accepted input and output nodes are:     
`LOSM, LO1, LO2, ..., LO12, \
NLOSM, NLO1, NLO2, ..., NLO12, \
cHHH0, cHHH2, cHHH5 \
LOSMfake2016, LO2fake2016, LO3fake2016, ..., LO13fake2016, \ 
LOSMfake2016WWgg, LO1fake2016WWgg, LO2fake2016WWgg, ..., LO12fake2016WWgg, \ 
LOSMfake2017, LO2fake2017, LO3fake2017, ..., LO12fake2017, \ 
LOSMfake2018, LO2fake2018, LO3fake2018, ..., LO12fake2018 `


## simple_example_usage.py     
This small script provides a simple usage example that can be easily  adapted to a specific analysis workflow.

## reweight_tree.py     
Is a reweighter tool that can be used in principle out-of-the-box. It reweights a set of
input TTrees, relative to a specific inputnode, to a set of outputnodes. The script contains
a set of extra options, including the possibility of production of validation plots. A usage example is:
` python reweight_tree.py \
--infilename /path/to/input.root \
--intreenames inputtreename1,inputtreename2,... \
--inputnode NLOcHHH5 \
--outputnodes NLOcHHH0,NLOSM \
--outdir /path/to/output/dir/ \
--mHHbranchname nameofgenMHHbranch \
--ValidationPlots \
--extrascaling 2. `       
