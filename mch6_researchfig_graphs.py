# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:07:52 2017

For making figs for research day - needs to have mch6-distAnalysis run first

@author: James Rig
"""

x = rats['mch6.5'].sessions['s15']

# Rename sipper events to distractors

distractors = x.output.Sir.onset

# Find distracted and non-distracted

lone_distractors = [val for idx, val in enumerate(distractors) if (val-distractors[idx-1] > 3)]

firstlick, distractedArray = jmf.distractedOrNot(lone_distractors, x.lickData['licks'])

# Use mastersnipper to make snips of distracted and non-distracted

snips_distracted = jmf.mastersnipper(x, lone_distractors)

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(3.2, 2.4))

jmfig.trialsMultShadedFig(ax[0], [x.siUVTrials[x.distractedArray],
                               x.siTrials[x.distractedArray]],
                          x.pps,
                          eventText = '')

jmfig.trialsMultShadedFig(ax[1], [x.siUVTrials[~x.distractedArray],
                               x.siTrials[~x.distractedArray]],
                          x.pps,
                          eventText = '')