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
lambda        = 10;
tol           = 1e-5;
maxiter       = 500;
dt            = 0.1;

[output,mask]=inpainting_harmonic(imagefilename,maskfilename,lambda,tol,maxiter,dt);
output_images(:,:,i)=output;
delete('image.png');delete('mask.png');

end

figure;imshow3D(images);
figure;imshow3D(output_images);