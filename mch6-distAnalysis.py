# -*- coding: utf-8 -*-
"""
Created on Wed May 24 11:09:00 2017

MCH6 analysis

@author: James Rig
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
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
    
    def __init__(self, data):      
        self.rat = data
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
        self.uvfile = datafolder + self.hrow['rat'] + 's' + self.hrow['session'] + '-v.mat'
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
        self.dataUV = self.output.resultUV
        
    def time2samples(self):
        tick = self.output.Tick.onset
        maxsamples = len(tick)*int(self.fs)
        if (len(self.data) - maxsamples) > 2*int(self.fs):
            print('Something may be wrong with conversion from time to samples')
            print(str(len(self.data) - maxsamples) + ' samples left over. This is more than double fs.')
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
        else:
            self.t2sMap = np.linspace(min(tick), max(tick), maxsamples)
            
    def event2sample(self, EOI):
        idx = (np.abs(self.t2sMap - EOI)).argmin()   
        return idx
    
    def check4events(self):        
        if hasattr(self.output, 'LiA'):
            self.licksA = self.output.LiA.onset
            self.offsetA = self.output.LiA.offset
        else:
            self.licksA = []
            self.offsetA = []
            
        if hasattr(self.output, 'LiB'):
            self.licksB = self.output.LiB.onset
            self.offsetB = self.output.LiB.offset
        else:
            self.licksB = []
            self.offsetB = []
            
    def sessionFig(self, ax):
        ax.plot(self.data, color='blue')
        try:
            ax.plot(self.dataUV, color='m')
        except:
            print('No UV data.')
        ax.set_xticks(np.multiply([0, 10, 20, 30, 40, 50, 60],60*self.fs))
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '60'])
        ax.set_xlabel('Time (min)')
        ax.set_title('Rat ' + self.rat + ': Session ' + self.session)

    def findlicks1bottle(self):
        if self.bottleB == 'empty':
            self.licks = self.output.LiA.onset
            self.offset = self.output.LiA.offset
            self.contents = self.bottleA
        elif self.bottleA == 'empty':
            self.licks = self.output.LiB.onset
            self.offset = self.output.LiB.offset
            self.contents = self.bottleB
        else:
            print('No empty bottle. Unable to assign licks. Please update metafile.')
           
def makeBehavFigs(x):
    # Initialize figure
    behavFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(left=0.10, right= 0.9, wspace=0.5, hspace = 0.7)
    plt.suptitle('Rat ' + x.rat + ': Session ' + x.session + ' (' + x.contents + ')')
    
    # Make figures of lick data from each bottle
    ax = plt.subplot(gs1[0, :])
    jmfig.sessionlicksFig(ax, x.lickData)        
    
    ax = plt.subplot(gs1[1, 0])
    jmfig.licklengthFig(ax, x.lickData)
    
    ax = plt.subplot(gs1[1, 1])
    jmfig.iliFig(ax, x.lickData)
    
    ax = plt.subplot(gs1[2, 0])
    jmfig.burstlengthFig(ax, x.lickData, color3rdbar=True)

    ax = plt.subplot(gs1[2, 1])
    jmfig.ibiFig(ax, x.lickData)
    
    #plt.tight_layout()
    pdf_pages.savefig(behavFig)

def makePhotoFigs(x):
    # Initialize photometry figure
    photoFig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    gs1 = gridspec.GridSpec(5, 2)
    gs1.update(left=0.125, right= 0.9, wspace=0.4, hspace = 0.8)
    plt.suptitle('Rat ' + x.rat + ': Session ' + x.session + ' (' + x.contents + ')')

    ax = plt.subplot(gs1[0, :])
    x.sessionFig(ax)

    ax = plt.subplot(gs1[1, 0])
    jmfig.trialsFig(ax, x.LiTrials, x.pps, eventText = 'First lick',
                    ylabel = 'Delta F / F0')
    
    ax = plt.subplot(gs1[1, 1])
    jmfig.trialsMultShadedFig(ax, [x.LiUVTrials, x.LiTrials], x.pps,
                              eventText = 'First lick')
    
    ax = plt.subplot(gs1[2, 0])
    jmfig.heatmapFig(ax, x.LiTrials, pps = x.pps)

    ax = plt.subplot(gs1[2, 1])
    data = jmf.nearestevents(x.lickData['rStart'], x.lickData['licks'])
    nLicks = [len(i) for i in data]
    sortOrder = np.argsort(nLicks)
    jmfig.heatmapFig(ax, x.LiTrials, pps = x.pps, sortEvs = nLicks)  
    
    if hasattr(x, 'distractedArray'):
        ax = plt.subplot(gs1[3, 0])
        jmfig.trialsMultShadedFig(ax, [x.siUVTrials[x.distractedArray],
                                       x.siTrials[x.distractedArray]],
                                  x.pps,
                                  eventText = '',
                                  title='Distracted trials (n=%i)' % np.sum(x.distractedArray))

        ax = plt.subplot(gs1[3, 1])
        jmfig.trialsMultShadedFig(ax, [x.siUVTrials[~x.distractedArray],
                                       x.siTrials[~x.distractedArray]],
                                  x.pps,
                                  eventText = '',
                                  title='Non-distracted trials')
        
        ax = plt.subplot(gs1[4, 0])
        jmfig.heatmapFig(ax, x.siTrials, pps = x.pps)
        
        ax = plt.subplot(gs1[4, 1])
        jmfig.heatmapFig(ax, x.siTrials, pps = x.pps, sortEvs = x.firstlick)  

#    plt.savefig(x.rat + '.eps', format='eps', dpi=1000)
    pdf_pages.savefig(photoFig)

# Read in metafile

metafileData, metafileHeader = jmf.metafilereader(userhome + '\\Documents\\GitHub\\MCH6\\mch6-dist-forMatPy.txt')
exptsuffix = '-dist'
includecol = 5

rats = {}

for i in metafileData:
    if int(i[includecol]) == 1:
        rowrat = str(i[2])
        if rowrat not in rats:
            rats[rowrat] = Rat(rowrat)
        rats[rowrat].loadsession(i, metafileHeader)

#for i in ['mch6.5']:
#    pdf_pages = PdfPages(i + '.pdf')
#    for j in ['s15']:
    
for i in rats:
    pdf_pages = PdfPages('R:\\DA_and_Reward\\jem64\\1705_MCH6\\output\\' + i + exptsuffix + '.pdf')
    for j in rats[i].sessions:
        
        print('Analysing rat ' + i + ' in session ' + j)
    
        # Load in data from .mat file (convert from Tank first using Matlab script)
        x = rats[i].sessions[j]
        x.loadmatfile()
        # Work out time to samples
        x.time2samples()       
        # Find out which bottles have TTLs/Licks associated with them
        x.check4events()
        x.findlicks1bottle()
        
        x.lickData = jmf.lickCalc(x.licks,
                                  offset = x.offset,
                                  burstThreshold = 0.50)
             
        bins = 200
    
        try:
            x.siTrials, x.pps = jmf.snipper(x.data, x.output.Sir.onset,
                                                t2sMap = x.t2sMap, fs = x.fs, bins=bins)
            x.siUVTrials, x.pps = jmf.snipper(x.dataUV, x.output.Sir.onset,
                                                t2sMap = x.t2sMap, fs = x.fs, bins=bins)
            x.firstlick, x.distractedArray = jmf.distractedOrNot(x.output.Sir.onset, x.lickData['licks'])
        except TypeError:
            print('Sir event is not present or not an array.')

        if len(x.lickData['rStart']) > 0:
            x.LiTrials, x.pps = jmf.snipper(x.data, x.lickData['rStart'], 
                                            t2sMap = x.t2sMap,
                                            fs = x.fs,
                                            bins=bins)
            x.LiUVTrials, x.pps = jmf.snipper(x.dataUV, x.lickData['rStart'], 
                                            t2sMap = x.t2sMap,
                                            fs = x.fs,
                                            bins=bins)        
        makeBehavFigs(x)       
        makePhotoFigs(x)

    pdf_pages.close()
"""
test = rats['mch6.6'].sessions['s4'].bottles['BottleB'].iliFig
rats['mch6.6'].sessions['s5'].hrow
http://python-guide-pt-br.readthedocs.io/en/latest/writing/style/

to add:
    histograms on session licks - 4
    individual points on interburstintervals - 5
    colours for different solutions - 3
    loop for behav figs - 2
    photometry data split for different events and coloured - 1
    heatmap for photometry
    lick checker (e.g. remove very very short ilis)
"""

