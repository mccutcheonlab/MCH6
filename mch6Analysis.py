# -*- coding: utf-8 -*-
"""
Created on Wed May 24 11:09:00 2017

MCH6 analysis

@author: James Rig
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import collections as col

import JM_general_functions as jmf
import JM_custom_figs as jmfig

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

plt.style.use('seaborn-muted')

import os

userhome = os.path.expanduser('~')
datafolder = userhome + '\\Dropbox\\Python\\matlab files\\'

class Rat(object):
    
    nRats = 0
    nSessions = 0
    
    def __init__(self, name):      
        self.rat = name
        self.sessions = {}
        
        Rat.nRats += 1
                
    def loadsession(self, data, header):
        self.session = 's'+str(data[4]) #should reference column of data with session number
        self.sessions[self.session] = Session(data, header, self.rat, self.session)
       
        Rat.nSessions += 1
        
class Session(object):
    
    def __init__(self, data, header, rat, session):
        self.hrow = {}
        for idx, col in enumerate(header):
            self.hrow[col] = data[idx]
        self.matlabfile = datafolder + self.hrow['rat'] + 's' + self.hrow['session'] + '.mat'
        self.medfile = datafolder + self.hrow['medfile']
        self.bottleA = self.hrow['bottleA']
        self.bottleB = self.hrow['bottleB']
        self.rat = str(rat)
        self.session = session
        
        self.bottles = {}
#        self.medfile = currentpath + 'MATLAB\\Experiments\\2015_casein\\behavfiles\\' + data[1]
#        self.lickdata = jmf.medfilereader(self.medfile, ['a', 'b', 'c', 'd'])
    def loadmatfile(self):
        a = sio.loadmat(self.matlabfile, squeeze_me=True, struct_as_record=False) 
        self.output = a['output']
        self.fs = self.output.fs
        self.data = self.output.result

    def loadUVfile(self):
        a = sio.loadmat(self.uvfile, squeeze_me=True, struct_as_record=False)
        b = a['output']
        self.dataUV = b.result
        
    def time2samples(self):
        tick = self.output.Tick.onset
        maxsamples = len(tick)*int(self.fs)
        if (len(self.data) - maxsamples) > 2*int(self.fs):
            print('Something may be wrong with conversion from time to samples')
            print(str(len(self.data) - maxsamples) + ' samples left over. This is more than double fs.')
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
        else:
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
            
#    def event2sample(self, EOI):
#        idx = (np.abs(self.t2sMap - EOI)).argmin()   
#        return idx
    
    def check4events(self):
        
        if hasattr(self.output, 'LiA'):
            self.ATrials = True
            self.licksA = self.output.LiA.onset
            self.offsetA = self.output.LiA.offset
        else:
            self.ATrials = False
            self.licksA = []
            self.offsetA = []
            
        if hasattr(self.output, 'LiB'):
            self.BTrials = True
            self.licksB = self.output.LiB.onset
            self.offsetB = self.output.LiB.offset
        else:
            self.BTrials = False
            self.licksB = []
            self.offsetB = []
            
        if hasattr(self.output, 'Sir'):
            self.sipper = self.output.Sir.onset
            self.sipper_off = self.output.Sir.offset
            
    def removephantomlicks(self):
        if self.ATrials == True:
            phlicks = jmf.findphantomlicks(self.licksA, self.sipper, postsipper=1.5)
            self.licksA = np.delete(self.licksA, phlicks)
            self.offsetA = np.delete(self.offsetA, phlicks)
    
        if self.BTrials == True:
            phlicks = jmf.findphantomlicks(self.licksB, self.sipper, postsipper=1.5)
            self.licksB = np.delete(self.licksB, phlicks)
            self.offsetB = np.delete(self.offsetB, phlicks)
            
    def sessionFig(self, ax):
        ax.plot(self.data)
        ylimits=ax.get_ylim()
        for x1,x2 in zip(self.sipper_off, self.sipper[1:]):
            ax.fill_between([x1*self.fs, x2*self.fs], ylimits[0], ylimits[1], color='grey', alpha=0.4)
        ax.set_xticks(np.multiply([0, 10, 20, 30, 40, 50, 60],60*self.fs))
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60'])
        ax.set_xlabel('Time (min)')
        ax.set_title('Rat ' + self.rat + ': Session ' + self.session)
        
    def sessionLicksFig(self, ax):
        if self.ATrials == True:
            ax.hist(self.licksA[1:], range(0, 3600, 60))
        if self.BTrials == True:
            ax.hist(self.licksB[1:], range(0, 3600, 60))
        ylimit = ax.get_ylim()
        ax.set_ylim([ylimit[0], ylimit[1]*1.4])
        yA = ylimit[1]*1.3
        yB = ylimit[1]*1.1
        if  self.ATrials == True:
            yraster = [yA for i in self.licksA[1:]]
            ax.scatter(self.licksA[1:], yraster, s=50, facecolors='none', edgecolors='grey')
        if self.BTrials == True:
            yraster = [yB for i in self.licksB[1:]]
            ax.scatter(self.licksB[1:], yraster, s=50, facecolors='none', edgecolors='grey')
       
        ax.set_xticks(np.multiply([0, 10, 20, 30, 40, 50, 60],60))
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60'])
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Licks per min')

        ylimits=ax.get_ylim()
        for x1,x2 in zip(self.sipper_off, self.sipper[1:]):
            ax.fill_between([x1, x2], ylimits[0], ylimits[1], color='grey', alpha=0.4)
        
        ax.set_xlabel('Time (s)')
        ax.set_title('Rat ' + self.rat + ': Session ' + self.session)
        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        xAB = ax.get_xlim()[1]
        ax.text(xAB, yA, self.bottleA, ha='right',va='center')
        ax.text(xAB, yB, self.bottleB, ha='right',va='center')
#        
    def makephotoTrials(self, bins, events, threshold=10):
        bgMAD = jmf.findnoise(self.data, self.randomevents,
                              t2sMap = self.t2sMap, fs = self.fs, bins=bins,
                              method='sum')          
        blueTrials, self.pps = jmf.snipper(self.data, events,
                                            t2sMap = self.t2sMap, fs = self.fs, bins=bins)        
    #    UVTrials, self.pps = jmf.snipper(self.dataUV, events,
    #                                        t2sMap = self.t2sMap, fs = self.fs, bins=bins)
        sigSum = [np.sum(abs(i)) for i in blueTrials]
        sigSD = [np.std(i) for i in blueTrials]
        noiseindex = [i > bgMAD*threshold for i in sigSum]
    
        return blueTrials, noiseindex
    #    return blueTrials, UVTrials, noiseindex
           
def makeBehavFigs(x):
    # Initialize figure
    behavFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(left=0.10, right= 0.9, wspace=0.5, hspace = 0.7)
    
    # Make figures of lick data from each bottle
    ax = plt.subplot(gs1[0, :])
    x.sessionLicksFig(ax)        
    
    
    if len(x.lickDataA['licks']) > 1:
        ax = plt.subplot(gs1[1, 0])
        jmfig.licklengthFig(ax, x.lickDataA, contents = x.bottleA)
        
        ax = plt.subplot(gs1[2, 0])
        jmfig.iliFig(ax, x.lickDataA, contents = x.bottleA)
        
        ax = plt.subplot(gs1[3, 0])
        jmfig.burstlengthFig(ax, x.lickDataA, contents = x.bottleA)

        ax = plt.subplot(gs1[4, 0])
        jmfig.ibiFig(ax, x.lickDataA, contents = x.bottleA)
    else:
        print('No licks for Bottle A in this session.')
        
    if len(x.lickDataB['licks']) > 1:
        ax = plt.subplot(gs1[1, 1])
        jmfig.licklengthFig(ax, x.lickDataB, contents = x.bottleB)
        
        ax = plt.subplot(gs1[2, 1])
        jmfig.iliFig(ax, x.lickDataB, contents = x.bottleB)
        
        ax = plt.subplot(gs1[3, 1])
        jmfig.burstlengthFig(ax, x.lickDataB, contents = x.bottleB)

        ax = plt.subplot(gs1[4, 1])
        jmfig.ibiFig(ax, x.lickDataB, contents = x.bottleB)        
    else:
        print('No licks for Bottle B in this session.')
    
    #plt.tight_layout()
    pdf_pages.savefig(behavFig)

def makePhotoFigs(x):
    # Initialize photometry figure
    photoFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)

    ax = plt.subplot(gs1[0, :])
    x.sessionFig(ax)

#    if x.siTrials != None:
    try:
        ax = plt.subplot(gs1[1, 0])
        jmfig.trialsFig(ax, x.siTrials, x.pps, eventText = 'Sipper Out')
        
        ax = plt.subplot(gs1[1, 1])
        jmfig.trialsShadedFig(ax, x.siTrials, x.pps, eventText = 'Sipper Out')        
    except:
        print('No no')
        
    if len(x.lickDataA['licks']) > 1: 
        ax = plt.subplot(gs1[2, 0])
        jmfig.trialsFig(ax, x.LiATrials, x.pps, eventText = x.bottleA)
        
        ax = plt.subplot(gs1[2, 1])
        jmfig.trialsShadedFig(ax, x.LiATrials, x.pps, eventText = x.bottleA)

    if len(x.lickDataB['licks']) > 1: 
        ax = plt.subplot(gs1[3, 0])
        jmfig.trialsFig(ax, x.LiBTrials, x.pps, eventText = x.bottleB)
        
        ax = plt.subplot(gs1[3, 1])
        jmfig.trialsShadedFig(ax, x.LiBTrials, x.pps, eventText = x.bottleB)
        
    if x.ATrials == True and x.BTrials == True:
        ax = plt.subplot(gs1[4, 0])
        jmfig.trialsMultShadedFig(ax, [x.LiATrials, x.LiBTrials], x.pps)
#
#    ax = plt.subplot(gs1[3, :])        
#    rats[i].sessions[j].trialsHeatMap(ax)
#    plt.savefig(userhome + '\\Dropbox\\Python\\figformasa.eps', format='eps', dpi=1000)
    pdf_pages.savefig(photoFig)

def makeHeatmapFigs(x):
    heatmapsFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(3, 2)
    gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
    
    ax = plt.subplot(gs1[0, 0])
    jmfig.heatmapFig(ax, x.siTrials, pps = x.pps)
    data = jmf.nearestevents(x.output.Sir.onset, x.lickDataA['licks'])
    jmfig.addevent2heatmap(ax, data, pps = x.pps)
    
    ax = plt.subplot(gs1[0, 1])
    data = jmf.nearestevents(x.output.Sir.onset, x.lickDataA['licks'])
    sortEv = jmf.findfirst(data)
    jmfig.heatmapFig(ax, x.siTrials, pps = x.pps, sortEvs = sortEv)

    jmfig.addevent2heatmap(ax, data, pps = x.pps)
    
    ax = plt.subplot(gs1[1, 0])
    jmfig.heatmapFig(ax, x.LiATrials, pps = x.pps)
    data = jmf.nearestevents(x.lickDataA['rStart'], x.output.Sir.onset) 
    jmfig.addevent2heatmap(ax, data, pps = x.pps)
    
    pdf_pages.savefig(heatmapsFig)

def checkphantomlicks(x):
    # This function was a debugging function designed to let me examine the lengths of 
    # phantom licks that occur as the sipper retracts from the box. It appears that
    # these licks are slightly shorter than most normal licks but not be enough, 
    # and not consitently enough for us to be able to easily remove them in the course
    # of lick clean up that we hope to perform as standard
    sipper = x.output.Sir.onset
    onset = x.licksA
    offset = x.offsetA
    
    length=[]
    for x in sipper:
        on = ([y for y in onset if y > x][0])
        off = ([y for y in offset if y > x][0])
        length.append(off-on)
    
    print(length)
    print(np.mean(length))
    
#    print(np.mean(x.lickDataA['licklength']))
#    print(np.std(x.lickDataA['licklength']))
#    
#    print(len(x.lickDataA['licklength']))
#    onelickind = [ind for ind, val in enumerate(x.lickDataA['licks'][:-1]) if
#                  (val - x.lickDataA['licks'][ind-1] > 0.25) &
#                  (x.lickDataA['licks'][ind+1] - val > 0.25)]
#    print(onelickind)
#
#    burstCount = col.Counter(x.lickDataA['bLicks'])
#    print(burstCount)
#    
#    onelicklength = x.lickDataA['licklength'][onelickind]
#    
    return onset, offset
    
    
# Read in metafile
metafileData, metafileHeader = jmf.metafilereader(userhome + '/Dropbox/Python/matlab files/mch6-forMatPy.txt')
exptsuffix = ''
includecol = 14

rats = {}

for i in metafileData:
    if int(i[includecol]) == 1:
        rowrat = str(i[2])
        if rowrat not in rats:
            rats[rowrat] = Rat(rowrat)
        rats[rowrat].loadsession(i, metafileHeader)

for i in ['mch6.6']:
    pdf_pages = PdfPages(userhome + '/Dropbox/Python/photometry/output-mch6/' + i + '.pdf')
    for j in ['s6']:
        
        print('Analysing rat ' + i + ' in session ' + j)
    
        x = rats[i].sessions[j]
        x.loadmatfile()

        x.time2samples()       

        x.check4events()
        x.lickDataA = jmf.lickCalc(x.licksA,
                           offset = x.offsetA)
        
        x.lickDataB = jmf.lickCalc(x.licksB,
                   offset = x.offsetB)
        
        x.removephantomlicks()
        
        bins=300
        x.randomevents = jmf.makerandomevents(120, max(x.output.Tick.onset)-120)
        x.bgTrials, x.pps = jmf.snipper(x.data, x.randomevents,
                                        t2sMap = x.t2sMap, fs = x.fs, bins=bins)
        
        x.siTrials, x.sinoise = x.makephotoTrials(bins, x.sipper_off)
        
        if x.ATrials == True:
            x.LiATrials, x.LiAnoise = x.makephotoTrials(bins, x.lickDataA['rStart'])
        
        if x.BTrials == True:
            x.LiBTrials, x.LiBnoise = x.makephotoTrials(bins, x.lickDataB['rStart'])

#        checkphantomlicks(x)
      
        makeBehavFigs(x)       
        makePhotoFigs(x)

#        makeHeatmapFigs(x)

    pdf_pages.close()
"""
http://python-guide-pt-br.readthedocs.io/en/latest/writing/style/

to add:
    histograms on session licks - 4
    individual points on interburstintervals - 5
    colours for different solutions - 3
    loop for behav figs - 2
    photometry data split for different events and coloured - 1
    heatmap for photometry
    lick checker (e.g. remove very very short ilis)
    check length of phantom licks
    
"""

