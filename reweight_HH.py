import os
import ROOT
from optparse import OptionParser
from array import array

def load_mHH_shapes(inputshapesLO, inputshapesNLO):
    mHH_shapes={}
    benchmarks = ["SM","1","2","3","4","5","6","7","8","9","10","11","12"]

    #load LO shapes
    infileLO = ROOT.TFile(inputshapesLO)
    for benchmark in benchmarks:
        for year in ["2016","2017"]:
            objectname = year+"_"+benchmark+"_LO_fake"
            benchmarkname = "LO"+benchmark
            if year=="2017":
                benchmarkname+="fake"
            mHH_shapes[benchmarkname] = infileLO.Get(objectname)
            mHH_shapes[benchmarkname].SetDirectory(0)
    infileLO.Close()

    #load NLO shapes
    infileNLO = ROOT.TFile(inputshapesNLO)
    for benchmark in benchmarks:
        objectname = "EFT_"+benchmark+"_NLO"
        benchmarkname = "NLO"+benchmark
        mHH_shapes[benchmarkname] = infileNLO.Get(objectname)
        mHH_shapes[benchmarkname].SetDirectory(0)
    #add also anomalous LO and NLO kl
    for extrabenchmark in ["cHHH0","cHHH2","cHHH5"]:
        for order in ["LO","NLO"]:
            objectname = "EFT_"+extrabenchmark+"_"+order
            benchmarkname = order+extrabenchmark
            mHH_shapes[benchmarkname] = infileNLO.Get(objectname)
            mHH_shapes[benchmarkname].SetDirectory(0)

    infileNLO.Close()

    return mHH_shapes

def ComputeSF(input_mHH, output_mHH):
    #################################################
    #################################################
    # TEMPFIX to preserve total number of events 
    # TO BE REMOVED with next version of histograms 
    input_mHH.Scale(1./input_mHH.Integral("width"))
    output_mHH.Scale(1./output_mHH.Integral("width"))
    #################################################
    #################################################
    mHH_SF=output_mHH.Clone()
    mHH_SF.Divide(input_mHH)
    mHH_SF.SetDirectory(0)
    return mHH_SF


parser = OptionParser()
parser.add_option("--infilename",      default="",                               help="Input file(s) for reweight")
parser.add_option("--intreenames",     default='*',                              help="Input tree(s) to be rereweighted")
parser.add_option("--mHHbranchname",   default="genMhh",                         help="Branch name containing the gen-mHH value")
parser.add_option("--inputnode",       default="",                               help="Node used for generation")
parser.add_option("--outputnodes",     default="NLOSM,NLO2,NLOcHHH2p45,NLOcHH5", help="Target nodes")
parser.add_option("--outdir",          default="./rew/",                         help="Output dir for plots")
parser.add_option("--inputshapesLO",   default="./reweight_HH_fake.root",        help="Input shapes for LO benchmarks")
parser.add_option("--inputshapesNLO",  default="./reweight_HH.root",             help="Input shapes for NLO benchmarks")
parser.add_option("--ValidationPlots", action="store_true", default=False,       help="Produce validation plots" )

(options,args)=parser.parse_args()

# load all the mHH benchmarks shapes for LO and NLO in mHH_shapes dictionary
# mHH_shapes should be used as mHH_shapes[<order><benchmark>]
# e.g. mHH_shapes["LOSM"], or mHH_shapes["LO6fake"] 
print "Loading mHH shapes"
mHH_shapes=load_mHH_shapes(options.inputshapesLO, options.inputshapesNLO) 

print "Opening ",options.infilename
infile = ROOT.TFile(options.infilename,"READ")
reduced_infilename = os.path.basename(options.infilename).replace(".root","")#to be used afterwards for outfile naming

intreenames = []
infolders = []
for fulltreename in options.intreenames.split(','):
    if '/' in fulltreename:
        infolders.append(fulltreename[:fulltreename.rfind('/')])
        intreenames.append(fulltreename[fulltreename.rfind('/')+1:])
    else:
        infolders.append("")
        intreenames.append(fulltreename)

for outputnode in options.outputnodes.split(','):
    print "reweight %s --> %s"%(options.inputnode,outputnode)

    # deriving SF for the requested inputnode-->outputnode 
    input_mHH = mHH_shapes[options.inputnode]
    output_mHH = mHH_shapes[outputnode]
    mHH_SF = ComputeSF(input_mHH,output_mHH)
    mHH_max = mHH_SF.GetXaxis().GetXmax()

    # opening outfile
    outfilename = "%s/%s_reweight_%s_to_%s.root"%(options.outdir,reduced_infilename,options.inputnode,outputnode)
    outfile = ROOT.TFile(outfilename,"RECREATE")
    
    #loop over trees
    for i_tree in range(0,len(intreenames)):
        #loading input tree
        treename=intreenames[i_tree]
        foldername=infolders[i_tree]
        intree = infile.Get("%s/%s"%(foldername,treename))
        mHH = array('f', [ 0 ])
        weight = array('f', [ 0 ])
        intree.SetBranchStatus("*",1)
        intree.SetBranchAddress(options.mHHbranchname,mHH)
        intree.SetBranchAddress("weight",weight)

        #prepare output tree
        outfile.cd()
        outdir=outfile.GetDirectory(foldername)
        if not outdir:
             outfile.mkdir(foldername)
        outfile.cd(foldername)
        outtree = intree.CloneTree(0)
        BMreweight = array('f', [ 0 ])
        outtree.Branch("BMreweight",BMreweight,'BMreweight/F')

        #loop over entries
        i=0;
        nEntries= intree.GetEntries()
        for i in range(0, nEntries):
            intree.GetEntry(i)
            if i%10000==0:
                print "reading entry ",i
            
            #compute reweight 
            reweight=0.
            if mHH[0]<250.:
                print "gen_mHH < 250 ... something is wrong"
            elif mHH[0]>=mHH_max:
                reweight = mHH_SF.GetBinContent( mHH_SF.GetXaxis().GetNbins() )#last bin
            else:
                reweight = mHH_SF.GetBinContent( mHH_SF.GetXaxis().FindBin(mHH[0]) )

            #update outfile weight
            weight[0] *= reweight
            BMreweight[0] = reweight
            #print mHH[0], weight[0], BMreweight[0]
            outtree.Fill()
        
        outtree.AutoSave()
    outfile.Close()
        
infile.Close()


#########################################################################
#########################################################################
# PLOTS
if options.ValidationPlots:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    infile = ROOT.TFile(options.infilename,"READ")
    for outputnode in options.outputnodes.split(','):
        print "Validation plots for %s --> %s"%(options.inputnode,outputnode)

        # retrieve input mHH, output mHH, and SF
        input_mHH = mHH_shapes[options.inputnode].Clone()
        output_mHH = mHH_shapes[outputnode].Clone()
        mHH_SF = ComputeSF(input_mHH,output_mHH)
        mHH_max = mHH_SF.GetXaxis().GetXmax()

        #draw options
        input_mHH.SetLineColor(1)
        input_mHH.SetLineWidth(2)
        output_mHH.SetLineColor(2)
        output_mHH.SetLineWidth(2)

        # re-open outfile
        outfilename = "%s/%s_reweight_%s_to_%s.root"%(options.outdir,reduced_infilename,options.inputnode,outputnode)
        outfile = ROOT.TFile(outfilename,"READ")
    
        #loop over trees
        c1 = ROOT.TCanvas()
        categorizedinputmHHs={}
        categorizedoutputmHHs={}
        categorizedinputmHHs["inclusive"] = ROOT.TH1F("inclusiveinputmHH","inclusiveinputmHH",100,0.,2000.)
        categorizedoutputmHHs["inclusive"] = ROOT.TH1F("inclusiveoutputmHH","inclusiveoutputmHH",100,0.,2000.)
        for i_tree in range(0,len(intreenames)):
            
            #loading input and output trees
            treename=intreenames[i_tree]
            foldername=infolders[i_tree]
            intree = infile.Get("%s/%s"%(foldername,treename))
            outtree = outfile.Get("%s/%s"%(foldername,treename))

            #build mHH histos
            intree.Draw("%s >> inputmHH_%s(100,0.,2000.)"%(options.mHHbranchname,treename),"(weight)","goff")
            categorizedinputmHHs[treename] = ROOT.gROOT.FindObject("inputmHH_%s"%treename)
            outtree.Draw("%s >> outputmHH_%s_%s(100,0,2000)"%(options.mHHbranchname,treename,outputnode),"(weight)","goff")
            categorizedoutputmHHs[treename] = ROOT.gROOT.FindObject("outputmHH_%s_%s"%(treename,outputnode))

            #fill inclusive distributions
            categorizedinputmHHs["inclusive"].Add(categorizedinputmHHs[treename])
            categorizedoutputmHHs["inclusive"].Add(categorizedoutputmHHs[treename])

        #draw mHH distributions
        for catname, categorizedinputmHH in categorizedinputmHHs.items():
            print "doing", catname
            categorizedoutputmHH = categorizedoutputmHHs[catname]

            #draw options
            categorizedinputmHH.SetLineColor(3)
            categorizedinputmHH.SetLineWidth(2)
            categorizedoutputmHH.SetLineColor(4)
            categorizedoutputmHH.SetLineWidth(2)

            # normalize mHH distributions expected for input and 
            # output benchmarks to the observed number in a given category  
            xmin = input_mHH.GetXaxis().GetXmin() 
            xmax = input_mHH.GetXaxis().GetXmax() 
            Nev_categorizedinput = categorizedinputmHH.Integral( categorizedinputmHH.GetXaxis().FindBin(xmin),
                                                                 categorizedinputmHH.GetXaxis().FindBin(xmax))
            input_mHH.Scale( Nev_categorizedinput / input_mHH.Integral() )
            Nev_categorizedoutput = categorizedoutputmHH.Integral( categorizedoutputmHH.GetXaxis().FindBin(xmin),
                                                                   categorizedoutputmHH.GetXaxis().FindBin(xmax))
            output_mHH.Scale( Nev_categorizedoutput / output_mHH.Integral() )
            
            #draw
            c1.cd()
            c1.Clear()
            categorizedinputmHH.Draw("hist")
            input_mHH.Draw("hist same")
            categorizedoutputmHH.Draw("hist same")
            output_mHH.Draw("hist same")
            categorizedinputmHH.GetXaxis().SetTitle("m_{HH}")

            leg = ROOT.TLegend(0.5,0.6,0.9,0.9)
            leg.AddEntry(categorizedinputmHH,"m_{HH} for selected ev. before rew.","l")
            leg.AddEntry(input_mHH,"m_{HH} for benchmark %s"%options.inputnode)
            leg.AddEntry(categorizedoutputmHH,"m_{HH} for selected ev. after rew.","l")
            leg.AddEntry(output_mHH,"m_{HH} for benchmark %s"%outputnode,"l")
            leg.Draw()
            title = ROOT.TLatex()
            title.SetNDC()
            title.DrawLatex(0.1,0.93,"#scale[0.6]{gen m_{HH} for cat. %s}"%catname)
            c1.Print("%s/mHH_from_%s_to_%s_category_%s.png"%(options.outdir,options.inputnode,outputnode,catname))
            c1.Print("%s/mHH_from_%s_to_%s_category_%s.pdf"%(options.outdir,options.inputnode,outputnode,catname))

        outfile.Close()

    infile.Close()

    
