clc;clear all;close all
load updated_dataset_1;
output_images=zeros(217,347,7);

for i=1:7
    

a=images(:,:,i);
bc=mask_images(:,:,i);
b=imcomplement(bc);
imwrite(a,'image.png')
imwrite(b,'mask.png')

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
output_images(:,:,i)=output;
delete('image.png');delete('mask.png');

end

figure;imshow3D(images);
figure;imshow3D(output_images);