import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
import numpy as np
import Util as utl

def plotVariableDistribution(histStore,varName,savePrefix=None,year='2018',lumi='58',lumiScale=False,sortByCounts=False):
    if savePrefix:
        f1=plt.figure(figsize=(8,8))
        ax1=plt.subplot(4,1,(1,3))
        ax0=plt.subplot(4,1,4)
        f2=plt.figure(figsize=(6,6))
        ax2=plt.subplot(4,1,(1,3))
        ax3=plt.subplot(4,1,4)
    else:
        f1=plt.figure(figsize=(14,6))
        ax1=plt.subplot(1,2,1)
        ax2=plt.subplot(4,2,(2,6))
        ax3=plt.subplot(4,2,8)
    
    hists=[]
    hists_integral=[]
    labels=[]
    histTot=None
    print("Drawing ",varName," for year ",year,"@ lumi = ",lumi)
    for ky  in histStore['bkg']:
        if ky not in utl.backgroundStackList:
            print(f"== skipping background : {ky} , since colurs are not defiend for the {ky} ! ")
            continue
    for ky  in utl.backgroundStackList: 
        if ky not in histStore['bkg']: 
            continue
        if not histTot:
            histTot=histStore['bkg'][ky].Clone()
        else:
            histTot.Add(histStore['bkg'][ky].Clone())
        
        hists.append(histStore['bkg'][ky].Clone())
        labels.append(ky)
        hists_integral.append(histStore['bkg'][ky].Integral())

        if lumiScale:
            hists[-1].Scale(float(lumi))

    if lumiScale:
        histTot.Scale(float(lumi))
    
    hists_toplot  = hists
    labels_toplot = labels
    if sortByCounts:
        ordering=np.argsort(hists_integral)
        hists_toplot=[hists[i] for i in ordering]
        labels_toplot=[labels_toplot[i] for i in ordering]
    else:
       labels_toplot=[]
       hists_toplot=[]
       for ky in utl.backgroundStackList:
            if ky in labels:
                idx=labels.index(ky)
                labels_toplot.append(labels[idx])
                hists_toplot.append(hists[idx])
    
    if len(labels)!=len(labels_toplot):
        print(f"{labels=} , {labels_toplot =} .. sorry they dont agree please update Util.backgroundStackList ! exiting ")
        exit(1)

    hep.histplot(hists_toplot,stack=True,label=labels_toplot,histtype='fill',ax=ax1)
    hep.histplot(histStore['data'],ax=ax1,label='Data',color='black',histtype='errorbar')
    ax1.legend(loc=0)
    ax1.semilogy()
    ax1.grid(alpha=0.1,color='k')
    ax1.set_xlabel(varName)
    
    rHistRaw=histStore['data'].Clone()
    rHistRaw.Divide(histTot)
    hep.histplot(rHistRaw,ax=ax0,color='magenta')
    #ax0.set_ylim([0.0,5])
    ax0.axhline(1.0,c='k')
    ax0.set_xlabel(varName)
    ax0.grid(alpha=0.1,color='k')

    histTotNorm=histTot.Clone() ;
    histStoreDataNorm=histStore['data'].Clone()

    if histTotNorm.Integral()==0.0:
        print("Total of bacgrounds giving Zero for variable ",varName,"  Exiting !")
        return
    else:
        histTotNorm.Scale(1.0/histTotNorm.Integral())
    if histTotNorm.Integral()==0.0:
        print("Total of Data giving Zero for variable ",varNamem," Exiting")
        return
    else:
        histStoreDataNorm.Scale(1.0/histStoreDataNorm.Integral())
    histTotNormSig=histStore['sig']['ggHHH'].Clone() ;histTotNormSig.Scale(1.0/histTotNormSig.Integral())

    hep.histplot(histTotNorm,ax=ax2,label='$\gamma(\gamma)$ + Jets',color='r')
    hep.histplot(histStoreDataNorm,ax=ax2,label='Data',color='b')
    hep.histplot(histTotNormSig,ax=ax2,label='ggHHH',color='darkorange')
    ax2.annotate("Shape Distribution",xy=(0.5, 0.5), xycoords='axes fraction')
    ax2.legend()
    ax2.set_xticklabels([])
    ax2.grid(alpha=0.1,color='k')

    rHist=histStoreDataNorm.Clone()
    rHist.Divide(histTotNorm)
    hep.histplot(rHist,ax=ax3,color='magenta')
    #ax3.set_ylim([0.0,5.0])
    ax3.axhline(1.0,c='k')
    ax3.set_xlabel(varName)
    ax3.grid(alpha=0.1,color='k')
    hep.cms.label("Preliminary",data=True,year=year,lumi=lumi,ax=ax1)
    hep.cms.label("Preliminary",data=True,year=year,lumi=lumi,ax=ax2)
    if savePrefix:
        f1.savefig(savePrefix+'/'+varName+'.png' ,bbox_inches='tight')
        f2.savefig(savePrefix+'/'+varName+'_shape.png' ,bbox_inches='tight')
    plt.close('all')
