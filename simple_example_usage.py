import reweight_HH 

inputnode="LO3fake2018" #string to define the input node
outputnodes=["NLOSM", "NLO1", "NLO2", "NLO3", "NLO4", "NLO5", "NLO6", "NLO7", "NLO8", "NLO9", "NLO10", "NLO11"] #collection of strings to define the output nodes
inputshapesLOfake="./shapes_v3/LOfake.root"#input file containing shapes for fake benchmarks
inputshapesLO="./shapes_v3/MadGraphLO.root"#input file containing shapes for LO benchmarks
inputshapesNLO="./shapes_v3/reweight_HH.root"#input file containing shapes for NLO benchmarks and NLO anomalous kl
extrascaling=1.# in this example I am considering only one inputnode. In case I want to reweight N input nodes (normalized to the same XS*lumi) to the same NLO benchmark, to then merge them, I have to scale each inputnode by an extra factor equal to N 

# define a reweighter object for each outputnode
HHreweighters = {}
for outputnode in outputnodes:
    HHreweighters[outputnode] = reweight_HH.HHreweighter(inputnode, outputnode, inputshapesNLO, inputshapesLO, inputshapesLOfake)
    
# simulate the input events
events=[]
events.append( {
    "weight":0.1,
    "mHH_gen":456.,
    "BMreweight":{},
} )
events.append( {
    "weight":0.08,
    "mHH_gen":568.,
    "BMreweight":{},
} )

for event in events:
    for outputnode in outputnodes:
        event["BMreweight"] [outputnode] = HHreweighters[outputnode].getWeight( event["mHH_gen"] ) / extrascaling
    
#done 
print events
