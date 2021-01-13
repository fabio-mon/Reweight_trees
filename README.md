Reweight samples.
=========
-Prepare flashgg output tree files.
-----------
The name of the sample should like this:
GluGluToHHTo2G2l2nu_node_1_2017.root

"{process}_{node_X}_{year}.root"

-Customize some variables inside the C programs.
------ 
![image](https://github.com/chuwang1/Reweight_trees/blob/main/Variables.png)
You need to provide: 
-Process name,
-year, 
-target_nodes(you can put more than one nodes at here.) 
-input_node(your sample's node).

Then you can change the vector of cats.
Set your categroy inside this cats vector.By default,there are two cats(HHWWggTag_2 and HHWWggTag_3).
I have added many systematics into the systematics vector.
If you have any other systematics, please add them into systematics vector,too.

Onece you complete these two steps.
then you can run this program.It will return a new root file and one validation plots.
The output file named {process}_{node_X}_{year}_reweighted.root









