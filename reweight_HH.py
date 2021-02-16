import os
import ROOT
#import copy
from optparse import OptionParser
from array import array

def load_mHH_shapes(inputshapesNLO, inputshapesLO, inputshapesLOfake):
    mHH_shapes={}
    benchmarks = ["SM","1","2","3","4","5","6","7","8","9","10","11","12"]

    # load LO FAKE shapes (WWgg 2016 shapes differ from bbgg!)
    # they do not need any scaling
    print "Loading LOfake shapes"
    infileFAKE = ROOT.TFile(inputshapesLOfake,"READ")
    for objkey in infileFAKE.GetListOfKeys():
        obj = objkey.ReadObj()
        if 'TH1' in obj.ClassName():
            mHH_shapes[obj.GetName()] = obj
            mHH_shapes[obj.GetName()].SetDirectory(0)
    infileFAKE.Close()

    # load LO shapes
    # they do not need any scaling
    print "Loading LO shapes"
    infile = ROOT.TFile(inputshapesLO,"READ")
    for benchmark in benchmarks:
        benchmarkname = "LO"+benchmark
        mHH_shapes[benchmarkname] = infile.Get(benchmarkname)
        mHH_shapes[benchmarkname].SetDirectory(0)
    infile.Close()

    # load NLO shapes
    # scaling by total XS extracted from xsections histogram
    print "Loading NLO shapes"
    infileNLO = ROOT.TFile(inputshapesNLO)
    h_xsections = infileNLO.Get("xsections")
    for benchmark in benchmarks:
        for order in ["NLO"]:
            objectname = "EFT_"+benchmark+"_"+order
            if benchmark=="8":
                objectname = "EFT_8a_"+order                
            benchmarkname = order+benchmark
            mHH_shapes[benchmarkname] = infileNLO.Get(objectname)
            mHH_shapes[benchmarkname].SetDirectory(0)
            XStot = h_xsections.GetBinContent( h_xsections.GetXaxis().FindBin(objectname) )
            if XStot==0.:
                raise RuntimeError("load_mHH_shapes : XStot for %s is found to be zero"%objectname)
            mHH_shapes[benchmarkname].Scale(1./XStot)
    #add also anomalous LO and NLO kl
    for extrabenchmark in ["cHHH0","cHHH2","cHHH5","SM"]:#I am overwriting NLO SM 
        for order in ["NLO"]:
            objectname = "EFT_"+extrabenchmark+"_"+order+"_V2"
            benchmarkname = order+extrabenchmark
            mHH_shapes[benchmarkname] = infileNLO.Get(objectname)
            mHH_shapes[benchmarkname].SetDirectory(0)
            #XStot = h_xsections.GetBinContent( h_xsections.GetXaxis().FindBin(objectname) )
            #if XStot==0.:
            #    raise RuntimeError("load_mHH_shapes : XStot for %s is found to be zero"%objectname)
            #mHH_shapes[benchmarkname].Scale(1./XStot)

    infileNLO.Close()

    return mHH_shapes

def ComputeSF(input_mHH, output_mHH):
    #find common mHH interval, and largest mHH interval
    commonmHHmin = max(input_mHH.GetXaxis().GetXmin(),output_mHH.GetXaxis().GetXmin())
    commonmHHmax = min(input_mHH.GetXaxis().GetXmax(),output_mHH.GetXaxis().GetXmax())
    mHHmin = min(input_mHH.GetXaxis().GetXmin(),output_mHH.GetXaxis().GetXmin())
    mHHmax = max(input_mHH.GetXaxis().GetXmax(),output_mHH.GetXaxis().GetXmax())
    print "Scale factors can be defined for %.2f < mHH < %.2f"%(commonmHHmin,commonmHHmax)
    print "Input shape integral for %.2f < mHH < %.2f is %.3f"%(commonmHHmin,commonmHHmax,input_mHH.Integral(input_mHH.GetXaxis().FindBin(commonmHHmin),input_mHH.GetXaxis().FindBin(commonmHHmax),"width"))
    print "Output shape integral for %.2f < mHH < %.2f is %.3f"%(commonmHHmin,commonmHHmax,output_mHH.Integral(output_mHH.GetXaxis().FindBin(commonmHHmin),output_mHH.GetXaxis().FindBin(commonmHHmax),"width"))
    stdbinwidth = 20.
    Nbin = int((mHHmax-mHHmin)/stdbinwidth)
    mHH_SF = ROOT.TH1F("SF_%s_%s"%(input_mHH.GetTitle(),output_mHH.GetTitle()),
                       "SF_%s_%s"%(input_mHH.GetTitle(),output_mHH.GetTitle()),
                       Nbin, mHHmin, mHHmax)

    mHH_SF.SetDirectory(0)
    for ibin in range(1,Nbin+1):
        mHH = mHH_SF.GetXaxis().GetBinCenter(ibin)

        if mHH<commonmHHmin:
            SF = output_mHH.GetBinContent(output_mHH.GetXaxis().FindBin(commonmHHmin)) / input_mHH.GetBinContent(input_mHH.GetXaxis().FindBin(commonmHHmin))
        elif mHH>=commonmHHmax:
            SF = output_mHH.GetBinContent(output_mHH.GetXaxis().FindBin(commonmHHmax-0.1)) / input_mHH.GetBinContent(input_mHH.GetXaxis().FindBin(commonmHHmax-0.1)) 
        else:
            output_mHH_value = output_mHH.GetBinContent(output_mHH.GetXaxis().FindBin(mHH))
            input_mHH_value = input_mHH.GetBinContent(input_mHH.GetXaxis().FindBin(mHH))
            if output_mHH_value<=0.:
                SF=0.
                if mHH>240.:
                    print "[WARNING]: output mHH shape is <= 0 for mHH=%.1f"%mHH
            elif input_mHH_value<=0.:
                SF=0.
                if mHH>240.:
                    print "[WARNING]: input mHH shape is <= 0 for mHH=%.1f"%mHH
            else:
                SF = output_mHH.GetBinContent(output_mHH.GetXaxis().FindBin(mHH)) / input_mHH.GetBinContent(input_mHH.GetXaxis().FindBin(mHH)) 

        mHH_SF.SetBinContent(ibin, SF)

    return mHH_SF

def InterpolateTail(h_mHH,xminfit):
    ROOT.gROOT.SetBatch(True)
    rp = ROOT.TFitResultPtr()
    fitfunc = ROOT.TF1("fitfunc_%s"%h_mHH.GetName(),"[0]*exp(-[1]*x)",xminfit,h_mHH.GetXaxis().GetXmax())
    fitfunc.SetParLimits(0,0.,100.)
    fitfunc.SetParLimits(1,0.,10.)
    fitfunc.SetParameters(10.,0.01)
    rp=h_mHH.Fit(fitfunc,"RS")
    fStatus = int(rp);
    if fStatus==4:
        print "[WARNING]: interpolation of %s did not converge"%h_mHH.GetName()
        return h_mHH
    c1 = ROOT.TCanvas()
    #h_mHH.Draw("E1")
    #c1.SaveAs("interpolation_%s.root"%h_mHH.GetName())
    
    normalization = h_mHH.Integral("width")
    ibin0 = h_mHH.GetXaxis().FindBin(xminfit)
    for ibin in range(ibin0,h_mHH.GetXaxis().GetNbins()+1):
        mHH = h_mHH.GetXaxis().GetBinCenter(ibin)
        h_mHH.SetBinContent(ibin,fitfunc.Eval(mHH))
    h_mHH.Scale( normalization/h_mHH.Integral("width") )

    return h_mHH
    

parser = OptionParser()
parser.add_option("--infilename",        default="",                               help="Input file(s) for reweight")
parser.add_option("--intreenames",       default='*',                              help="Input tree(s) to be rereweighted")
parser.add_option("--mHHbranchname",     default="genMhh",                         help="Branch name containing the gen-mHH value")
parser.add_option("--inputnode",         default="",                               help="Node used for generation")
parser.add_option("--outputnodes",       default="NLOSM,NLO2,NLOcHHH2p45,NLOcHH5", help="Target nodes")
parser.add_option("--outdir",            default="./rew/",                         help="Output dir for plots")
parser.add_option("--inputshapesLOfake", default="./shapes_v3/LOfake.root",       help="Input shapes for fake benchmarks")
parser.add_option("--inputshapesLO",     default="./shapes_v3/MadGraphLO.root",   help="Input shapes for LO benchmarks")
parser.add_option("--inputshapesNLO",    default="./shapes_v3/reweight_HH.root",  help="Input shapes for NLO benchmarks")
parser.add_option("--extrascaling",      default=1.,   type="float",               help="extra scaling factor to include for reweight")
parser.add_option("--ValidationPlots",   default=False, action="store_true",       help="Produce validation plots" )
parser.add_option("--InterpolateInputTail",   default=False, action="store_true",       help="Smooth mHH tail of the input through interpolation" )
parser.add_option("--xmininputinterp",      default=500.,   type="float",               help="xmin for input mHH interpolation")
parser.add_option("--InterpolateOutputTail",   default=False, action="store_true",       help="Smooth mHH tail of the output through interpolation" )
parser.add_option("--xminoutputinterp",      default=500.,   type="float",               help="xmin for output mHH interpolation")
(options,args)=parser.parse_args()

# load all the mHH benchmarks shapes for LO and NLO in mHH_shapes dictionary
# mHH_shapes should be used as mHH_shapes[<order><benchmark>]
# e.g. mHH_shapes["LOSM"], or mHH_shapes["LO6fake"] 
print "Loading mHH shapes"
mHH_shapes=load_mHH_shapes(options.inputshapesNLO, options.inputshapesLO, options.inputshapesLOfake)

print "Opening ",options.infilename
infile = ROOT.TFile(options.infilename,"READ")
reduced_infilename = os.path.basename(options.infilename).replace(".root","")#to be used afterwards for outfile naming

intreenames = []
infolders = []
for fulltreename in options.intreenames.split(','):
    if '/' in fulltreename:
        infolders.append(fulltreename[:fulltreename.rfind('/')+1])
        intreenames.append(fulltreename[fulltreename.rfind('/')+1:])
    else:
        infolders.append("")
        intreenames.append(fulltreename)

for outputnode in options.outputnodes.split(','):
    print "reweight %s --> %s"%(options.inputnode,outputnode)

    # deriving SF for the requested inputnode-->outputnode 
    input_mHH = mHH_shapes[options.inputnode].Clone()
    output_mHH = mHH_shapes[outputnode].Clone()
    if options.InterpolateInputTail:
        input_mHH = InterpolateTail(input_mHH,options.xmininputinterp)
    if options.InterpolateOutputTail:
        output_mHH = InterpolateTail(output_mHH,options.xminoutputinterp)
    mHH_SF = ComputeSF(input_mHH,output_mHH)
    mHH_min = mHH_SF.GetXaxis().GetXmin()
    mHH_max = mHH_SF.GetXaxis().GetXmax()

    # opening outfile
    outfilename = "%s/%s_reweight_%s_to_%s.root"%(options.outdir,reduced_infilename,options.inputnode,outputnode)
    outfile = ROOT.TFile(outfilename,"RECREATE")
    
    #loop over trees
    for i_tree in range(0,len(intreenames)):

        #loading input tree
        treename=intreenames[i_tree]
        foldername=infolders[i_tree]
        intree = infile.Get("%s%s"%(foldername,treename))
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
            #if i%10000==0:
            #    print "reading entry ",i
            
            #compute reweight 
            reweight=1.
            if mHH[0]<250.:
                print "gen_mHH < 250 ... something is wrong"
            elif mHH[0]<=mHH_min:
                reweight = mHH_SF.GetBinContent(1)#first bin
            elif mHH[0]>=mHH_max:
                if(mHH[0]>2000): #over 2000 GeV do not reweight to avoid crazy values
                    reweight=1.
                else:
                    reweight = mHH_SF.GetBinContent( mHH_SF.GetXaxis().GetNbins() )#last bin
            else:
                reweight = mHH_SF.GetBinContent( mHH_SF.GetXaxis().FindBin(mHH[0]) )

            #update outfile weight
            weight[0] *= (reweight/options.extrascaling)
            BMreweight[0] = (reweight/options.extrascaling)
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
        if options.InterpolateInputTail:
            input_mHH = InterpolateTail(input_mHH,options.xmininputinterp)
        if options.InterpolateOutputTail:
            output_mHH = InterpolateTail(output_mHH,options.xminoutputinterp)
        mHH_SF = ComputeSF(input_mHH,output_mHH)
        commonmHHmin = max(input_mHH.GetXaxis().GetXmin(),output_mHH.GetXaxis().GetXmin())
        commonmHHmax = min(input_mHH.GetXaxis().GetXmax(),output_mHH.GetXaxis().GetXmax())

        #draw options
        input_mHH.SetLineColor(15)
        input_mHH.SetLineWidth(1)
        output_mHH.SetLineColor(860+8)
        output_mHH.SetLineWidth(1)

        # re-open outfile
        outfilename = "%s/%s_reweight_%s_to_%s.root"%(options.outdir,reduced_infilename,options.inputnode,outputnode)
        outfile = ROOT.TFile(outfilename,"READ")
    
        #loop over trees
        c1 = ROOT.TCanvas("c","c",500,500)
        categorizedinputmHHs={}
        categorizedoutputmHHs={}
        categorizedinputmHHs["inclusive"] = ROOT.TH1F("inclusiveinputmHH","inclusiveinputmHH",100,0.,2000.)
        categorizedoutputmHHs["inclusive"] = ROOT.TH1F("inclusiveoutputmHH","inclusiveoutputmHH",100,0.,2000.)
        for i_tree in range(0,len(intreenames)):
            
            #loading input and output trees
            treename=intreenames[i_tree]
            foldername=infolders[i_tree]
            intree = infile.Get("%s%s"%(foldername,treename))
            outtree = outfile.Get("%s%s"%(foldername,treename))

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

            if not "inclusive" in catname: continue

            print "doing", catname
            categorizedoutputmHH = categorizedoutputmHHs[catname]

            #draw options
            categorizedinputmHH.SetMarkerColor(1)
            categorizedinputmHH.SetMarkerStyle(20)
            categorizedinputmHH.SetMarkerSize(0.5)
            #categorizedinputmHH.SetLineWidth(2)
            categorizedoutputmHH.SetMarkerColor(4)
            categorizedoutputmHH.SetMarkerStyle(21)
            categorizedoutputmHH.SetMarkerSize(0.4)
            #categorizedoutputmHH.SetLineWidth(2)

            # normalize mHH distributions expected for input and 
            # output benchmarks to the observed number in a given category  
            Nev_categorizedinput = categorizedinputmHH.Integral( categorizedinputmHH.GetXaxis().FindBin(input_mHH.GetXaxis().GetXmin()+0.1),
                                                                 categorizedinputmHH.GetXaxis().FindBin(input_mHH.GetXaxis().GetXmax()-0.1))
            input_mHH.Scale( Nev_categorizedinput / input_mHH.Integral() )
            Nev_categorizedoutput = categorizedoutputmHH.Integral( categorizedoutputmHH.GetXaxis().FindBin(output_mHH.GetXaxis().GetXmin()+0.1),
                                                                   categorizedoutputmHH.GetXaxis().FindBin(output_mHH.GetXaxis().GetXmax()-0.1))
            output_mHH.Scale( Nev_categorizedoutput / output_mHH.Integral() )

            #ratioplot
            input_mHH_ratio = categorizedinputmHH.Clone()
            output_mHH_ratio = categorizedoutputmHH.Clone()
            for ibin in range(1,input_mHH_ratio.GetXaxis().GetNbins()):
                mHH=input_mHH_ratio.GetXaxis().GetBinCenter(ibin)
                input_mHH_value=input_mHH.GetBinContent(input_mHH.GetXaxis().FindBin(mHH))
                if input_mHH_value==0:
                    input_mHH_ratio.SetBinContent(ibin,0.)
                else:
                    obs_value = input_mHH_ratio.GetBinContent(ibin)
                    obs_value_err = input_mHH_ratio.GetBinError(ibin)
                    input_mHH_ratio.SetBinContent(ibin,obs_value/input_mHH_value)
                    input_mHH_ratio.SetBinError(ibin,obs_value_err/input_mHH_value)
                output_mHH_value=output_mHH.GetBinContent(output_mHH.GetXaxis().FindBin(mHH))
                if output_mHH_value==0:
                    output_mHH_ratio.SetBinContent(ibin,0.)
                else:
                    obs_value = output_mHH_ratio.GetBinContent(ibin)
                    obs_value_err = output_mHH_ratio.GetBinError(ibin)
                    output_mHH_ratio.SetBinContent(ibin,obs_value/output_mHH_value)
                    output_mHH_ratio.SetBinError(ibin,obs_value_err/output_mHH_value)

            #draw
            c1.cd()
            c1.Clear()
            pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1.0);
            pad1.SetBottomMargin(0.03);
            pad1.Draw();               
            pad1.cd(); 

            categorizedinputmHH.Draw("P")
            input_mHH.Draw("hist same")
            categorizedoutputmHH.Draw("P same")
            output_mHH.Draw("hist same")

            categorizedinputmHH.GetYaxis().SetRangeUser(0.,1.15*max(categorizedinputmHH.GetMaximum(),
                                                                    categorizedoutputmHH.GetMaximum()))
            categorizedinputmHH.GetXaxis().SetLabelSize(0.);
            categorizedinputmHH.GetXaxis().SetTitleSize(0.);


            leg = ROOT.TLegend(0.5,0.6,0.9,0.9)
            leg.AddEntry(categorizedinputmHH,"m_{HH} for selected ev. before rew.","p")
            leg.AddEntry(input_mHH,"m_{HH} for benchmark %s"%options.inputnode,"l")
            leg.AddEntry(categorizedoutputmHH,"m_{HH} for selected ev. after rew.","p")
            leg.AddEntry(output_mHH,"m_{HH} for benchmark %s"%outputnode,"l")
            leg.Draw()
            title = ROOT.TLatex()
            title.SetNDC()
            title.DrawLatex(0.1,0.93,"#scale[1.]{%s from %s to %s}"%(catname,options.inputnode,outputnode))
            title.DrawLatex(0.6,0.3,"#scale[0.8]{#frac{N_{ev}(after)}{N_{ev}(before)} = %.3f}"%(
                categorizedoutputmHH.Integral(0,-1)/categorizedinputmHH.Integral(0,-1)))
                #categorizedoutputmHH.Integral()/categorizedinputmHH.Integral()))

            c1.cd()
            pad2 = ROOT.TPad("pad2","pad2", 0,0.,1,0.3);
            pad2.SetTopMargin(0.03);
            pad2.SetBottomMargin(0.4);
            pad2.SetGridy();
            pad2.Draw();
            pad2.cd();
            input_mHH_ratio.Draw("P")
            output_mHH_ratio.Draw("P SAME")
            input_mHH_ratio.GetXaxis().SetTitle("m_{HH}")
            input_mHH_ratio.GetYaxis().SetRangeUser(0.,2.)
            input_mHH_ratio.GetYaxis().SetNdivisions(505);
            input_mHH_ratio.GetXaxis().SetTitleFont(43);
            input_mHH_ratio.GetYaxis().SetTitleFont(43);
            input_mHH_ratio.GetXaxis().SetTitleSize(25);
            input_mHH_ratio.GetYaxis().SetTitleSize(22);
            input_mHH_ratio.GetXaxis().SetTitleOffset(4.);
            input_mHH_ratio.GetYaxis().SetTitleOffset(1.85);
            input_mHH_ratio.GetXaxis().SetLabelFont(43);
            input_mHH_ratio.GetYaxis().SetLabelFont(43);
            input_mHH_ratio.GetXaxis().SetLabelSize(17);
            input_mHH_ratio.GetYaxis().SetLabelSize(17);



            c1.Print("%s/mHH_from_%s_to_%s_category_%s.png"%(options.outdir,options.inputnode,outputnode,catname))
            c1.Print("%s/mHH_from_%s_to_%s_category_%s.pdf"%(options.outdir,options.inputnode,outputnode,catname))

        outfile.Close()

    infile.Close()

    
