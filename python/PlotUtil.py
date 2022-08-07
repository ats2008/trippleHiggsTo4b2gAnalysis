import matplotlib.pyplot as plt
import ROOT 
import numpy as np
import os,copy
import PlotUtil as putil
import itertools as itrTools

# kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
# kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
# kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

import array
def getDefaultCanvas():
    
    defaultCanvasParams=dict({"CMS":{},"TAG":{},"EnLumi":{}})
    defaultCanvasParams["CMS"]["text"]         = "CMS"
    defaultCanvasParams["CMS"]["textSize"]     = 0.05
    defaultCanvasParams["CMS"]["textFont"]     = 61   
    defaultCanvasParams["CMS"]["xpos"]         = 0.15   
    defaultCanvasParams["CMS"]["ypos"]         = 0.98+0.004
    
    # for the "preliminary"
    defaultCanvasParams["TAG"]["text"]         = "internal"
    defaultCanvasParams["TAG"]["textSize"]     = 0.76*defaultCanvasParams["CMS"]["textSize"]  
    defaultCanvasParams["TAG"]["textFont"]     = 52   
    defaultCanvasParams["TAG"]["xpos"]         = 0.15  + 0.12
    defaultCanvasParams["TAG"]["ypos"]         = 0.98  - 0.004
    
    # for the "2018A [ 1 fb^{-1}] 13 TeV" label
    defaultCanvasParams["EnLumi"]["text"]       = "               13 TeV"
    defaultCanvasParams["EnLumi"]["textSize"]  = 0.76*defaultCanvasParams["CMS"]["textSize"]  
    defaultCanvasParams["EnLumi"]["textFont"]  = 42   
    defaultCanvasParams["EnLumi"]["xpos"]      = 0.85   
    defaultCanvasParams["EnLumi"]["ypos"]      = 0.95

    return defaultCanvasParams

def getCanvas2018A():
    cPars=getDefaultCanvas()
 #   cPars["EnLumi"]["text"]="   13 TeV"
    return cPars

def getCanvasParams(tag=None):
    cPars=getDefaultCanvas()
    print("Canvas Tag : ", tag )
    if tag=="2018Full":
        cPars["EnLumi"]["xpos"]       =  0.961
        cPars["EnLumi"]["text"]       = "SingleEG 2018 (58.9fb^{-1}) 13 TeV"
        return cPars
    if tag=="Run3MC":
        cPars["EnLumi"]["xpos"]       =  0.961
        cPars["EnLumi"]["text"]       = "Run 3 MC, 14 TeV"
        return cPars
    if tag=="GEN":
        cPars["EnLumi"]["xpos"]       =  0.961
        cPars["EnLumi"]["text"]       = "GEN , 13 TeV"
        return cPars
    if tag=="Run3_0p900":
        cPars["EnLumi"]["xpos"]       =  0.961
        cPars["EnLumi"]["text"]       = "Run 3,900 GeV"
        return cPars
    if tag=="Run2_2018":
        cPars["EnLumi"]["xpos"]       =  0.961
        cPars["EnLumi"]["text"]       = "2018, 59 fb"
        return cPars
    if tag=="Run2MC2018":
        cPars["EnLumi"]["xpos"]       =  0.961
        cPars["EnLumi"]["text"]       = "2018, MC"
        return cPars
    if tag=="Run2MC2017":
        cPars["EnLumi"]["xpos"]       =  0.961
        cPars["EnLumi"]["text"]       = "2017, MC"
        return cPars
    if tag!=None:
        cPars["EnLumi"]["text"]       = tag
    return cPars

def getFitParams(hist0):
    mean  = hist0.GetMean()
    stdDev= hist0.GetStdDev()
    minX=mean-4*stdDev
    maxX=mean+4*stdDev
    hist0.Fit("gaus","Q0","",minX,maxX);
    f=hist0.GetFunction("gaus");
    return [f.GetParameter(1),f.GetParameter(2)]

class Plot:
    def __init__(self, **args):
        self.name        = args.get("Name", "plot")
        self.legend      = args.get("Legend",None)
        self.histo       = args.get("Histo", None)
        self.fit         = args.get("Fit", None)
        self.markerColor = args.get("MarkerColor", ROOT.kBlue+2)
        self.markerStyle = args.get("MarkerStyle", 22)
        self.markerSize = args.get("MarkerSize", 5)
        self.lineColor   = args.get("LineColor", ROOT.kBlue+2)
        self.lineStyle   = args.get("LineStyle", 1)
        self.lineWidth   = args.get("LineWidth", 1)
        self.options   = args.get("Options", "")
        self.drawLegend   = args.get("DrawLegend",False)
        self.drawLine     = args.get("DrawLine", False)
        self.histo.SetName(self.name+"_histo")
        self.normalize =args.get("Normalize",False)
        self.scaleFactor =args.get("ScaleFactor",None)
        self.doFit =args.get("doFit",False)

class PlotCanvas:
    def __init__(self, **args):
        self.name  = ""
        self.plots = []
        self.rates = []
        self.plotDir =  args.get("prefix","plots/")
        self.yRange = (0.0, 1.0)
        self.xRange = (0, 100)
        self.yRange2 = (0.0, 1.0)
        self.xTitle = args.get("xTitle","E_{T}^{e offl} [GeV]")
        self.yTitle = args.get("yTitle","a.u")
        self.yTitle2 = args.get("yTitle","a.u")
        self.legendPosition = (0.4,0.2,0.9,0.6)
        self.descPosition = (0.6,0.58)
        self.desc = args.get("desc",["L1 Efficiency"])
        self.logx = args.get("logx",False)
        self.logy = args.get("logy",False)
        self.logy2 = args.get("logy2",False)
        self.canvasParams = args.get("canvasPars",getDefaultCanvas())    
        self.setPlotStyle()

    def addRate(self, rate):
        self.rates.append(rate)
    def addPlot(self, plot):
        self.plots.append(plot)
    def clearPlots(self):
        self.plots=[]


    def plot(self):
        xmax=0.0
        xmin=1e9
        xmax_=1e9
        xmin_=0.0
        ymax=0.0
        ymin=1e9
        ymax_=1e9
        ymin_=0.0
        qunats=np.zeros(2)
        probs=np.zeros(2)
        probs[0]=0.90
        probs[1]=0.95
        for plot in self.plots:
            histo = plot.histo
#             print("Scale Factor  : before  = ",plot.scaleFactor)    
            if plot.normalize:
                if plot.scaleFactor==None:
                    plot.scaleFactor=1.0/histo.Integral()
                else:
                    plot.scaleFactor/=histo.Integral()
#                 print("Integral = = ",1.0/histo.Integral())    
            if plot.scaleFactor!=None:
                histo.Scale(plot.scaleFactor)
#             print("Scaled the histogram by a factor of : ",plot.scaleFactor)
#             print("             Integral after Scaling : ",histo.Integral())
            
            if self.yRange[1]=='auto':
                ymax_=histo.GetMaximum()
                if ymax < ymax_*1.1:
                    ymax=ymax_*1.1
            if self.yRange[0]=='auto': 
                ymin_=histo.GetMinimum()
                if ymin > ymin_*0.8:
                    ymin=ymin_*0.8
            
            if self.xRange[1]=='auto':
                xmax_=histo.GetXaxis().GetBinCenter( histo.GetNbinsX() ) + 0.5*histo.GetBinWidth(histo.GetNbinsX())
#                 xmax_=hist.GetQuantiles(2,probs,qunats)
#                 xmax_=qunats[1]*1.2
#                 print(xmax_,xmax)
                if abs(xmax) < abs(xmax_):
                    xmax=xmax_
            if self.xRange[0]=='auto': 
                xmin_=histo.GetXaxis().GetBinCenter(1) - 0.5*histo.GetBinWidth(1)
                print(xmin_,xmin)
                if abs(xmin) > abs(xmin_):
                    xmin=xmin_
        
        if self.yRange[0]=='auto':
            self.yRange=(ymin,self.yRange[1])
        if self.yRange[1]=='auto':
            self.yRange=(self.yRange[0],ymax)
        if self.xRange[0]=='auto':
            self.xRange=(xmin,self.xRange[1])
        if self.xRange[1]=='auto':
            self.xRange=(self.xRange[0],xmax)

        print("Setting  X-RANGE  : ",self.xRange)
        print("Setting  Y-RANGE  : ",self.yRange)
        
#         print("xmax xmin : ",xmin,xmax)
        canvas = ROOT.TCanvas("c_"+self.name, self.name, 700, 800)
        canvas.SetGrid()
        hDummy = ROOT.TH1F("hDummy_"+self.name, self.name, 1, self.xRange[0],self.xRange[1])

        #print("ymin_,ymax_",ymin_,ymax_,"ymin , ymax - > " ,ymin,ymax)
            #print( "\tY : Min - > Max  : " ,histo.GetMinimum()," -> ", histo.GetMaximum())
        hDummy.SetAxisRange(self.yRange[0], self.yRange[1], "Y")
        hDummy.SetXTitle(self.xTitle)
        hDummy.SetYTitle(self.yTitle)
        hDummy.Draw("C")
        if self.logx : canvas.SetLogx()
        if self.logy : canvas.SetLogy()
        self.logz= True
        if self.logz : canvas.SetLogz()
        
        
        CMSbox=self.getCMSBox()
        extraTextBox=self.getExtraTextBox()
        lumibox=self.getLumiBox()


        selbox=[]
        n=0
        for inx in range(len(self.desc)):
            selbox.append(ROOT.TLatex  (self.descPosition[0], self.descPosition[1] -n*0.04, self.desc[n]))
            selbox[-1].SetNDC()
            selbox[-1].SetTextSize(0.035)
            selbox[-1].SetTextFont(12)
            selbox[-1].SetTextAlign(13)
            n+=1
        for plot in self.plots:
            if(plot.doFit):
                histo = plot.histo
                fitP  = getFitParams(histo)
                strFit='('+'{0:0.3f}'.format(fitP[0])+ ' , '+'{0:0.3f}'.format(fitP[1]) +')'
                selbox.append(ROOT.TLatex  (self.descPosition[0], self.descPosition[1] -n*0.04, strFit))
                selbox[-1].SetNDC()
                selbox[-1].SetTextSize(0.035)
                selbox[-1].SetTextFont(12)
                selbox[-1].SetTextColor(plot.lineColor)
                selbox[-1].SetTextAlign(13)
                n+=1

        #Line legend
        legend = ROOT.TLegend(self.legendPosition[0],self.legendPosition[1],self.legendPosition[2],self.legendPosition[3])
        legend.SetTextFont(42)
        legend.SetFillColor(0)

        for plot in self.plots:
            histo = plot.histo
            print("Ploting hist : ",histo.GetName())
            print( "\tX : Min - > Max  : " ,histo.GetBinCenter(1)," -> ", histo.GetBinCenter(histo.GetNbinsX())," [ ",histo.GetBinCenter(2)-histo.GetBinCenter(1)," ]" )
            print( "\tY : Min - > Max  : " ,histo.GetMinimum()," -> ", histo.GetMaximum())
            
            print(plot.markerSize,plot.markerStyle)
            histo.SetMarkerStyle(plot.markerStyle)
            histo.SetMarkerSize(plot.markerSize)
            histo.SetMarkerColor(plot.markerColor)
            histo.SetLineStyle(plot.lineStyle)
            histo.SetLineWidth(plot.lineWidth)
            histo.SetLineColor(plot.lineColor)
            #histo.Draw("same e")
#             print(plot.options)
            histo.Draw("same "+plot.options)
#             print( histo.GetOption())

            if(plot.legend):
                legend.AddEntry(histo, plot.legend, "pe")
                legend.SetTextSize(0.025) 
            if(plot.drawLegend) :  legend.Draw()
        
        scale = (self.yRange[1] -self.yRange[0])/(self.yRange2[1] - self.yRange2[0])
        y20=self.yRange2[0]

        if(self.logy2):
            scale = (self.yRange[1] -self.yRange[0])/(ROOT.log(self.yRange2[1]) - ROOT.log(self.yRange2[0]) )
            y20 = ROOT.log(self.yRange2[0])


        #####   begin :            For plotting a graph in the Second Axis         #####

        if len(self.rates) >0 :
            axis=None
            canvas.SetRightMargin(0.18);
            if(self.logy2):
                axis = ROOT.TGaxis(self.xRange[1],self.yRange[0],self.xRange[1], self.yRange[1],self.yRange2[0],self.yRange2[1],510,"+LG")
            else :
                axis = ROOT.TGaxis(self.xRange[1],self.yRange[0],self.xRange[1], self.yRange[1],self.yRange2[0],self.yRange2[1],510,"+L")
            axis.SetLineColor(ROOT.kRed);
            axis.SetTextColor(ROOT.kRed);
            axis.SetTitle(self.yTitle2)
            axis.SetTitleColor(ROOT.kBlack)
            axis.SetTitleOffset(1.2)
            axis.Draw();

        for rate in self.rates:
            histo = rate.histo
            if self.logy2:
                for i in range(histo.GetNbinsX()):
                    content=histo.GetBinContent(i+1)
                    scaledContent = ( ROOT.log(content + 1e-20 )  - y20 )* scale
                    histo.SetBinContent(i+1,scaledContent)
                    #print("\t ",content," -> ",scaledContent,"  = ( ",ROOT.log(content + 1e-20),"  - ",y20," )* ",scale)
            else:
                for i in range(histo.GetNbinsX()):
                    content=histo.GetBinContent(i+1)
                    scaledContent = (content  - y20 )* scale
                    histo.SetBinContent(i+1,scaledContent)
            
            histo.SetMarkerStyle(rate.markerStyle)
            histo.SetMarkerColor(rate.markerColor)
            histo.SetLineColor(rate.markerColor)
            histo.Draw("pl same "+ rate.options)
            legend.AddEntry(histo, rate.legend, "pe")
            legend.SetTextSize(0.025) 
            if(rate.drawLegend) :  legend.Draw()
        
        #####   end :            For plotting a graph in the Second Axis         #####
            
        CMSbox.Draw()
        extraTextBox.Draw()
        lumibox.Draw()
        for selb in selbox:
            selb.Draw()
        
#         canvas.Print(self.plotDir+"/"+self.name+".pdf", "pdf")
        canvas.Print(self.plotDir+"/"+self.name+".png", "png")

        return canvas

    def getCMSBox(self,pad=None):
        if pad!=None:
            pad.cd()
        par=self.canvasParams["CMS"]
        CMSbox= ROOT.TLatex  (par["xpos"],par["ypos"], par["text"])
        CMSbox.SetNDC()
        CMSbox.SetTextSize(par["textSize"])
        CMSbox.SetTextFont(par["textFont"])
        CMSbox.SetTextColor(ROOT.kBlack)
        CMSbox.SetTextAlign(13) 
        return CMSbox
   
    def getExtraTextBox(self):
        par=self.canvasParams["TAG"]
        extraTextBox = ROOT.TLatex  (par["xpos"], par["ypos"], par["text"])
        extraTextBox.SetNDC()
        extraTextBox.SetTextSize(par["textSize"])
        extraTextBox.SetTextFont(par["textFont"])
        extraTextBox.SetTextColor(ROOT.kBlack)
        extraTextBox.SetTextAlign(13)
        
        return extraTextBox

    def getLumiBox(self):
        par=self.canvasParams["EnLumi"]
        lumibox = ROOT.TLatex  ( par["xpos"], par["ypos"], par["text"])
        lumibox.SetNDC()
        lumibox.SetTextAlign(31)
        lumibox.SetTextSize(par["textSize"])
        lumibox.SetTextFont(par["textFont"])
        lumibox.SetTextColor(ROOT.kBlack)
        
        return lumibox
        
    def setPlotStyle(self):
        ROOT.gROOT.SetStyle("Plain")
        ROOT.gStyle.SetOptStat()
        ROOT.gStyle.SetOptFit(0)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetFrameLineWidth(2)
        ROOT.gStyle.SetPadBottomMargin(0.13)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadTopMargin(0.06)
        ROOT.gStyle.SetPadRightMargin(0.05)

        ROOT.gStyle.SetLabelFont(42,"X")
        ROOT.gStyle.SetLabelFont(42,"Y")
        ROOT.gStyle.SetLabelSize(0.04,"X")
        ROOT.gStyle.SetLabelSize(0.04,"Y")
        ROOT.gStyle.SetLabelOffset(0.01,"Y")
        ROOT.gStyle.SetTickLength(0.02,"X")
        ROOT.gStyle.SetTickLength(0.02,"Y")
        ROOT.gStyle.SetLineWidth(2)
        ROOT.gStyle.SetTickLength(0.02 ,"Z")

        ROOT.gStyle.SetTitleSize(0.1)
        ROOT.gStyle.SetTitleFont(42,"X")
        ROOT.gStyle.SetTitleFont(42,"Y")
        ROOT.gStyle.SetTitleSize(0.05,"X")
        ROOT.gStyle.SetTitleSize(0.05,"Y")
        ROOT.gStyle.SetTitleOffset(1.1,"X")
        ROOT.gStyle.SetTitleOffset(1.4,"Y")
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPaintTextFormat("3.2f")
        ROOT.gROOT.ForceStyle()

col=[ 0,1,3,4,6,ROOT.kOrange,ROOT.kMagenta+2,9,11,12,13,14,15  ]

def getDefaultPlotParams(col=1,marker=22,markerSize=1):
    return {'Legend':"Data : 2018, Unpacked",
                          'MarkerColor' : col,
                          'LineColor':col,
                          'LineWidth':2,
                          'LineStyle':9,
                          'MarkerStyle':marker,
                          'MarkerSize':markerSize,
                          'DrawLegend':True,
                          'Option':'pc'
            } 

def getDefaultPlot(prefix='plots/',name='defaultName',cPars=getCanvasParams('GEN')):
#     plot=putil.PlotCanvas(prefix=prefix,canvasPars=cPars)
    plot=PlotCanvas(prefix=prefix,canvasPars=cPars)
    
    plot.name   = name
    plot.xRange = (3.0,50.0)  
    plot.yRange = (0.0,1.02) 
    plot.yTitle = "#epsilon"  
    plot.xTitle = "E_{T}^{offl.}"
    plot.desc   = ["L1 EG Isolation"] 
    plot.legendPosition = (0.6,0.5,0.90,0.65)
    plot.descPosition   = (0.6,0.75)     
    plot.logx = False
    plot.logy = False
    return plot
