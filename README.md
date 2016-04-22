ngsParalog
==========

ngsParalog is a program for detecting cryptic copy number variation from population-level, next generation sequencing (NGS) data. ngsParalog implements a likelihood method for estimating the probability of mismapping reads due to paralogs while modeling aspects of NGS data, which makes it effective at even very low sequencing depths. The program reads in pileup format data and outputs per site likelihood ratios of duplication. These likelihood ratios are asymptotically distributed as a 50/50 mixture of a Chi-Square(1 d.f.) and a point mass at zero.

###Installation

Download ngsParalog with the command:

% git clone https://github.com/tplinderoth/ngsParalog

To install ngsParalog run the commands:

% cd ngsParalog
% make

To remove compilation products, run:

% cd ngsParalog
% make clean

###Author Details

Written by Tyler Linderoth
Contact: tylerp.linderoth@gmail.com 
