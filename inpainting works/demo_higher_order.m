clc;clear all;close all
inpainted_output=zeros(217,347,7);
load updated_dataset_1.mat
for i=1:7
    
lambda=mask(:,:,i);%this is the mask.....%%% for this algorithm this mask should be black
g=images(:,:,i);%this is the inpainted data
[m,n]=size(g);
g=double(g);
g = g./max(max(g));
clims=[0 1];
%Parameters
h1=1;
h2=1;
T=100;
dt=0.1;
lambda0=10;
lambda = lambda.*lambda0;
% figure(1); imagesc(g, clims); colormap(gray)

% laplacian = L1*u+u*L2
L1=(1/(h1^2))*(diag(-2*ones(m,1)) + diag(ones(m-1,1),1) + diag(ones(m-1,1),-1));
L1(1,1)=-1/h1^2; L1(m,m)=-1/h1^2;
L2=(1/(h2^2))*(diag(-2*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
L2(1,1)=-1/h2^2; L2(n,n)=-1/h2^2;

%intialization of u:
u=g;
% Explicit time stepping scheme for the heat equation
for t=1:T
    u=u+dt*(L1*u+u*L2+lambda.*(g-u));
%     figure(2);
%     imagesc(u, clims), axis equal; axis off; colormap(gray)
%     title(['Inpainted image as solution to heat equation after '...
%         num2str(t) ' timesteps']);
%     pause(0.01)
end

inpainted_output(:,:,i)=u;


end
figure;imshow3D(inpainted_output);title('output data');
figure;imshow3D(images);title('input data with artificial lesions');