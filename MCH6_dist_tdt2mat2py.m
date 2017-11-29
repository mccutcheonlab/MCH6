%% Convert MCH6 files from TDT to Python via Matlab

clear all; close all;

folder = 'C:\Users\James Rig\Dropbox\Python\mch6-dist\';
metafile = strcat(folder,'mch6-dist-forMatPy.txt');

fid = fopen(metafile);
C = textscan(fid, '%s %s %s %d %f %d %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

dataFolder = 'R:\DA_and_Reward\kp259\MCH6\MCH6 tanks\';
savefolder = 'C:\Users\James Rig\Dropbox\Python\matlab files\'

for i = 1:size(C{1,1},1)
    if C{1,6}(i) == 1 % checks to see if it is to be included or not
        TDTfile = char(strcat(dataFolder,C{1,1}(i)));
        rat = char(C{1,3}(i));
        session = strcat('s',num2str(C{1,5}(i)));
        tdt2mat2py(TDTfile,rat,session,0,savefolder)
    end
end

