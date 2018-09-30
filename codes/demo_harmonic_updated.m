clc;clear all;close all
load updated_dataset_408;
output_images=zeros(274,334,5);
for i=1:5
    
a=total_lesioned_image(:,:,180+i);
bc=total_binary_mask(:,:,180+i);

%b=imcomplement(bc);
imwrite(a,'image.png')
imwrite(bc,'mask.png')

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
figure;imshow3D(total_lesioned_image(:,:,181:185));
figure;imshow3D(output_images);
figure;imshow3D(L);title('our algorithm')
