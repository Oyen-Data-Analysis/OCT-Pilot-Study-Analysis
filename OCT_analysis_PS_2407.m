%% Volume Fraction, Surface Area per Volume, Pore Size analysis of OCT images for pilot study

clear;close all;clc

% Enter the full path and file name of binary image to analyze 
fileName = '/Users/adriennescott/Library/CloudStorage/Box-Box/Pilot Study Data/RR_Analysis/FGR003-1/segmented/3_P1_B2_NP_3D_230316_1_149-399segmented.tif'; 


%% Imports Data
fprintf('Currently reading file %s\n',fileName)
temptiff = tiffreadVolume(fileName);
imgStack = squeeze(temptiff(:,:,:,1));

imgStack = imgStack >= 1;

% %smoothing 
N = round(size(imgStack,1)/50);
kernel = ones(N, N) / N^2; % Create kernel for smoothing (was ^2)
blurryImage = convn(double(imgStack), kernel, 'same'); 
final_mask_image = blurryImage > 0.5;
A =bwareaopen(final_mask_image,50); % generate smooth segmented image

final_mask_image = final_mask_image(:,:,1:250);

%% Quantification of SA/V, PoreSizes, and Volume Fraction

quant=regionprops3(final_mask_image,"SurfaceArea","Extent","Volume"); 

total_vol_pixels=sum(quant.Volume); %in pixels

s=size(final_mask_image);

volume_fraction= total_vol_pixels/(s(1)*s(2)*s(3))

total_surface_area=sum(quant.SurfaceArea); %pixelssquared

SA_V_pixels=total_surface_area/total_vol_pixels;

SA_V_mm_1=(SA_V_pixels/((4*4*2.09)^(1/3)))*1000 %1/mm 

[pore_size] = Main_pore_f(final_mask_image, fileName);

avg_pore_size_um=mean(pore_size)

function [pore_size]=Main_pore_f(bw,FileName)

bw2 = bw(1:200,1:200,1:150);
A=bw2;

vol_ps(A)
L=segment_ps(A);
surf_ps(L)
Res=4; % spatial resolution 
REG=regionprops(L);

Equiv_Rad=([REG(:).Area].*3./4./pi).^(1/3)*Res;
name=erase(FileName,'.tif');
savename1=strcat(name,'_Pores.png');
saveas(gcf,savename1)

% plotting the pore size distribution
figure; histogram(Equiv_Rad,20,'Normalization','Probability'); 
xlabel('Pore radius (micron)');
ylabel('Probability')
title('Pore size distribution')

%plotting the fitted pore size distribution
figure; histfit(Equiv_Rad,20,'lognormal'); 
xlabel('Pore radius (micron)');
ylabel('Probability')
title('Pore size distribution')

pore_size=Equiv_Rad;


end

function vol_ps(A)
figure;
S=size(A);
A=1-A;
B=ones(S+2);
B(2:end-1,2:end-1,2:end-1)=A;
isosurface(B); axis equal tight;
end

function surf_ps(A)
A=double(A);
S=size(A);
S=permute(S,[2,1,3]);
figure;
[X,Y,Z]=meshgrid(1:S(1),1:S(2),1:S(3));
xslice=[1,S(1)];
yslice=[1,S(2)];
zslice=[1,S(3)];
h=slice(X,Y,Z,A,xslice,yslice,zslice); axis equal tight
for I=1:6
    h(I).EdgeColor='none';
end
end

function [L]=segment_ps(A)
A=double(A);
A=1-A;
B=bwdist(1-A);
B=imgaussfilt3(B,2);
L=watershed(-B);
L=double(L).*double(A);
end


