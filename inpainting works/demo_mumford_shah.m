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
maxiter       = 100; 
tol           = 1e-7;
param.lambda  = 10^5; % weight on data fidelity (should usually be large).
param.alpha   = 100;  % regularisation parameters \alpha.
param.beta    = 100;  % regularisation parameters \beta.
param.epsilon = 0.02; % accuracy of Ambrosio-Tortorelli approximation of the edge set.

[output,mask]=inpainting_mumford_shah(imagefilename,maskfilename,maxiter,tol,param);
output_images(:,:,i)=output;
delete('image.png');delete('mask.png');

end

figure;imshow3D(images);
figure;imshow3D(output_images);