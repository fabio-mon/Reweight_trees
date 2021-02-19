############################################################################### 
# class to perform reweight from any input to any output HH node (LO or NLO)
# It consists of a 1D reweight of the mHH variable
#
# see simple_example_usage.py for a simple usage example that can be easily 
# adapted to a specific analysis workflow
#    
# instead see reweight_tree.py for a out-of-the-box reweight of a given 
# input TTree to an output TTree
############################################################################### 

import os
import ROOT
#import copy
from optparse import OptionParser
from array import array

class HHreweighter:

    def __init__(self,inputnode, outputnode, inputshapesNLO, inputshapesLO, inputshapesLOfake):
        self.stdbinwidth = 20.
        self.inputnode = inputnode
        self.outputnode = outputnode
        self.mHH_shapes = self.load_mHH_shapes(inputshapesNLO, inputshapesLO, inputshapesLOfake)
        self.mHH_SF = self.ComputeSF()
        self.mHH_min = self.mHH_SF.GetXaxis().GetXmin()
        self.mHH_max = self.mHH_SF.GetXaxis().GetXmax()
        
    def inputshape(self):
        return self.mHH_shapes[self.inputnode].Clone()

    def outputshape(self):
        return self.mHH_shapes[self.outputnode].Clone()

    def PrintAvailableShapes(self):
        print "Available input/output nodes are"
        print self.mHH_shapes.keys()

    def getWeight(self,mHH):
        if mHH<250.:
            print "gen_mHH < 250 ... something is wrong"
        elif mHH<=self.mHH_min:
            reweight = self.mHH_SF.GetBinContent(1)#first bin
        elif mHH>=self.mHH_max:
            if(mHH>2000): #over 2000 GeV do not reweight to avoid crazy values
                reweight=1.
            else:
                reweight = self.mHH_SF.GetBinContent( self.mHH_SF.GetXaxis().GetNbins() )#last bin
        else:
            reweight = self.mHH_SF.GetBinContent( self.mHH_SF.GetXaxis().FindBin(mHH) )
        return reweight


    def load_mHH_shapes(self, inputshapesNLO, inputshapesLO, inputshapesLOfake):
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
            objectname = "EFT_"+benchmark+"_NLO"
            if benchmark=="8":
                objectname = "EFT_8a_NLO"
            benchmarkname = "NLO"+benchmark
            mHH_shapes[benchmarkname] = infileNLO.Get(objectname)
            mHH_shapes[benchmarkname].SetDirectory(0)
            XStot = h_xsections.GetBinContent( h_xsections.GetXaxis().FindBin(objectname) )
            if XStot==0.:
                raise RuntimeError("load_mHH_shapes : XStot for %s is found to be zero"%objectname)
            mHH_shapes[benchmarkname].Scale(1./XStot)
        #add also anomalous LO and NLO kl
        for extrabenchmark in ["cHHH0","cHHH2","cHHH5","SM"]:#I am overwriting NLO SM 
            objectname = "EFT_"+extrabenchmark+"_NLO_V2"
            benchmarkname = "NLO"+extrabenchmark
            mHH_shapes[benchmarkname] = infileNLO.Get(objectname)
            mHH_shapes[benchmarkname].SetDirectory(0)

        infileNLO.Close()

        return mHH_shapes


    def ComputeSF(self):
        input_mHH=self.mHH_shapes[self.inputnode]
        output_mHH=self.mHH_shapes[self.outputnode]

        #find common mHH interval, and largest mHH interval
        commonmHHmin = max(input_mHH.GetXaxis().GetXmin(),output_mHH.GetXaxis().GetXmin())
        commonmHHmax = min(input_mHH.GetXaxis().GetXmax(),output_mHH.GetXaxis().GetXmax())
        mHHmin = min(input_mHH.GetXaxis().GetXmin(),output_mHH.GetXaxis().GetXmin())
        mHHmax = max(input_mHH.GetXaxis().GetXmax(),output_mHH.GetXaxis().GetXmax())
        print "Scale factors can be defined for %.2f < mHH < %.2f"%(commonmHHmin,commonmHHmax)
        print "Input shape integral for %.2f < mHH < %.2f is %.3f"%(commonmHHmin,commonmHHmax,input_mHH.Integral(input_mHH.GetXaxis().FindBin(commonmHHmin),input_mHH.GetXaxis().FindBin(commonmHHmax),"width"))
        print "Output shape integral for %.2f < mHH < %.2f is %.3f"%(commonmHHmin,commonmHHmax,output_mHH.Integral(output_mHH.GetXaxis().FindBin(commonmHHmin),output_mHH.GetXaxis().FindBin(commonmHHmax),"width"))
        Nbin = int((mHHmax-mHHmin)/self.stdbinwidth)
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

