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

snips = jmf.mastersnipper(x, lone_distractors, preTrial=3, trialLength=10)

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(5.8, 3.5))

ax[0].axis('off')
ax[1].axis('off')

jmfig.trialsMultShadedFig(ax[0], [snips['uv'][distractedArray],
                               snips['blue'][distractedArray]],
                          x.pps,
                          eventText = '',
                          linecolor=['xkcd:silver', 'blue'], errorcolor=['xkcd:light grey', 'xkcd:light blue'])

jmfig.trialsMultShadedFig(ax[1], [snips['uv'][~distractedArray],
                               snips['blue'][~distractedArray]],
                          x.pps,
                          eventText = '',
                          linecolor=['xkcd:silver', 'xkcd:charcoal'], errorcolor=['xkcd:light grey', 'xkcd:grey'])

ax[0].plot([0,0], [0.02, 0.04], c='k')
ax[0].text(-10, 0.03, '2% \u0394F', verticalalignment='center', horizontalalignment='right')

plt.savefig('R:\\DA_and_Reward\\jem64\\1705_MCH6\\output\\photo-dist.eps')