clc;clear all;close all

load result_p5_t1.mat
[a,b]=size(mask);
rot_output_data=zeros(a,b);
output_images=zeros(a,b);
p=zeros(a,b);

rot_output_bin=zeros(a,b);
output_images_bin=zeros(a,b);

a=lesioned_image;
bc=mask;
imwrite(a,'image.png')
imwrite(bc,'mask.png')

imagefilename = 'image.png';
maskfilename  = 'mask.png';
% PARAMETERS
tol           = 1e-5;
maxiter       = 50;
dt            = 0.1;
param.M       = 40; % number of steps of the inpainting procedure;
param.N       = 2;  % number of steps of the anisotropic diffusion;
param.eps     = 1e-10;

% accuracy of Ambrosio-Tortorelli approximation of the edge set.

[output,mask]=inpainting_transport(imagefilename,maskfilename,maxiter,tol,dt,param);
output_images=output;




m.i=original_image;
rot_output_data=m.i ;
rot_output_data=mat2gray(rot_output_data);
th=0.5;

rot_output_bin=im2bw(rot_output_data,th);
   output_images_bin=im2bw(output_images,th);
   
rot_output_data=mat2gray(rot_output_data);
[Jaccard,Dice,rfp,rfn]=sevaluate(rot_output_bin,output_images_bin)


delete('image.png','mask.png','masked_harmonic.png','output_harmonic.png','log_harmonic.log');