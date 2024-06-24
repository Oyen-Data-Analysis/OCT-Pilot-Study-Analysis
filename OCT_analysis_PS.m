%% Volume Fraction, Surface Area per Volume, Pore Size analysis of OCT images for pilot study

clear;close all;clc

% Enter the full path and file name of binary image to analyze 
datadir = 'path/filename'; 


%% Imports Data
fprintf('Currently reading file %s\n',currentfilename)
temptiff = tiffreadVolume(fileName);
imgStack = squeeze(temptiff(:,:,:,1));

imgStack = imgStack >= 1;

%% Data Calculation
quant=regionprops3(imgStack,"SurfaceArea","Extent","Volume"); 

quant.Properties.VariableNames

