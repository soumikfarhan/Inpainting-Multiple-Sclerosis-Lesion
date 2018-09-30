% Copyright (C) 2011-2016 Simone Parisotto, Carola-Bibiane Schoenlieb
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  Authors:
%  Simone Parisotto (email: sp751 at cam dot ac dot uk)
%  Carola-Bibiane Schoenlieb (email: cbs31 at cam dot ac dot uk)
%      
%  Address:
%  Cambridge Image Analysis
%  Centre for Mathematical Sciences
%  Wilberforce Road
%  Cambridge CB3 0WA
%  United Kingdom
%  
%  Date:
%  February, 2016
%%

function [u_end,u_start]=inpainting_transport(imagefilename,maskfilename,maxiter,tol,dt,param)

%% ----------------- CREATE A log FILE WHERE TO STORE RESULTS IN txt FORMAT
logfilename = 'log_transport.log';
if exist(logfilename,'file')
    delete(logfilename);
end
fileID = fopen(logfilename,'w');

%% ------------------------------------ IMPORT THE CLEAN INPUT AND THE MASK
iminfo = imfinfo(imagefilename);
input  = im2double(imread(imagefilename));
% check if grayscale/truecolor dimension of image grey/colour
colors = size(input,3);

mask   = im2double(imread(maskfilename));
mask = double(mat2gray(mask)==0); % characteristic function for the intact part of the image
if size(mask,3)==1 && colors>1
    mask = repmat(mask,[1,1,colors]);
end

%% ---------------------------------------------------- INITIALIZATION OF u
u_start = (mask).*input + ~mask;
u       = reshape(u_start,iminfo.Height*iminfo.Width,colors);

%% ---------------------------------------------- GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;

%% --------------------------------------------------------------- GRADIENT
% FORWARD AND BACKWARD
d1i_forward  = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[0,1],iminfo.Height,iminfo.Height)/h1;
d1j_forward  = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[0,1],iminfo.Width,iminfo.Width)/h2;
d1i_backward = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[-1,0],iminfo.Height,iminfo.Height)/h1;
d1j_backward = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[-1,0],iminfo.Width,iminfo.Width)/h2;

% FOR BOUNDARY CONDITIONS
d1i_forward(end,:) = 0;
d1j_forward(end,:) = 0;
d1i_backward(1,:)  = 0;
d1j_backward(1,:)  = 0;

matrices.Dif  = kron(speye(iminfo.Width),d1i_forward);
matrices.Djf  = kron(d1j_forward,speye(iminfo.Height));
matrices.Dib  = kron(speye(iminfo.Width),d1i_backward);
matrices.Djb  = kron(d1j_backward,speye(iminfo.Height));

% CENTERED
d1i_centered  = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[-1,1],iminfo.Height,iminfo.Height)/(2*h1);
d1j_centered  = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[-1,1],iminfo.Width,iminfo.Width)/(2*h2);
% border not taken into account
d1i_centered([1 end],:) = 0;
d1j_centered([1 end],:) = 0;
matrices.Dic  = kron(speye(iminfo.Width),d1i_centered);
matrices.Djc  = kron(d1j_centered,speye(iminfo.Height));

%% -------------------------------------------------------------- LAPLACIAN
d2i = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Height))/h1^2;
d2j = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Width))/h2^2;
% PERIODIC BOUNDARY CONDITIONS
d2i(1,end) = 1/h1^2;
d2i(end,1) = 1/h1^2;
d2j(end,1) = 1/h2^2;
d2j(1,end) = 1/h2^2;
matrices.L = (kron(speye(iminfo.Width),d2i)+kron(d2j,speye(iminfo.Height)));

% - LAPLACIAN second derivatives with other bc
d2i = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Height))/h1^2;
d2j = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Width))/h2^2;
d2i([1 end],:) = 0;
d2j([1 end],:) = 0;
matrices.M2i = kron(speye(iminfo.Width),d2i);
matrices.M2j = kron(d2j,speye(iminfo.Height));
matrices.M2ij = matrices.Djc*matrices.Dic;

%% ------------------------------------------------------------ FREE MEMORY
clear d1i_forward d1j_forward
clear d1i_backward d1j_backward
clear d1i_centered d1j_centered
clear d2i d2j

%% DIFFUSION CONSTANT $g_{epsilon}$ within the small epsilon-strip around the inpainting domain:
SE           = strel('ball',6,0,0);
maskeps      = 1-mask;
lambdaeps    = reshape(imdilate(maskeps,SE)-maskeps,iminfo.Height*iminfo.Width,colors);
channel_mask = reshape(mask,iminfo.Height*iminfo.Width,colors);
maskeps      = reshape(maskeps,iminfo.Height*iminfo.Width,colors);
geps         = maskeps;
for k=1:colors
    for t=1:5
        % just interpolate g_{epsilon} with a few steps of linear diffusion within the strip
        geps(:,k) = geps(:,k) + dt*(matrices.L*geps(:,k)) + lambdaeps(:,k).*(maskeps(:,k)-geps(:,k));
    end
end

%% -------------------------------------------------------------- ALGORITHM

% INITIALIZATION
Dlap   = zeros(iminfo.Height*iminfo.Width,2);
normal = zeros(iminfo.Height*iminfo.Width,2);

% FOR EACH COLOR CHANNEL
for k=1:colors
    
    % ANISOTROPIC DIFFUSION PREPROCESSING STEP
    u(:,k) = anisodiff(u(:,k),dt,param.eps,geps(:,k),1,matrices);
    
    % ITERATION
    for iter = 1:maxiter
        % param.M steps of the inpainting procedure:
        for m=1:param.M
            lap         = matrices.L*u(:,k);
            Dlap(:,1)   = matrices.Dic*lap;
            Dlap(:,2)   = matrices.Djc*lap;
            normal(:,1) = matrices.Dic*u(:,k);
            normal(:,2) = matrices.Djc*u(:,k);
            
            normal  = bsxfun(@rdivide,normal,sqrt(sum(normal.^2,2) + param.eps));
            
            beta    = Dlap(:,1).*(-normal(:,2)) + Dlap(:,2).*normal(:,1) ;
            betapos = (beta>0);
            
            uxf = matrices.Dif*u(:,k);
            uxb = matrices.Dib*u(:,k);
            uyf = matrices.Djf*u(:,k);
            uyb = matrices.Djb*u(:,k);
            
            slopelim = betapos .* sqrt(min(uxb,0).^2 + max(uxf,0).^2 + min(uyb,0).^2 + max(uyf,0).^2)...
                + (~betapos) .* sqrt(max(uxb,0).^2 + min(uxf,0).^2 + max(uyb,0).^2 + min(uyf,0).^2);
            
            update = beta .* slopelim;
            % OPTIONAL:
            % Nonlinear scaling of the equation proposed by Bertalmio in 
            % his thesis. It might cause instabilities when choosing dt 
            % too large.
            %         signo = sign(update);
            %         update = signo.*sqrt(sqrt(signo.*update));
            
            % UPADTE ONLY PIXELS INSIDE THE INPAINTING DOMAIN (BY MASK)
            u(:,k) = u(:,k) + dt * ~channel_mask(:,k).*update;
        end
        
        % param.N steps of anisotropic diffusion
        un = anisodiff(u(:,k),dt,param.eps,geps(:,k),param.N,matrices);
        
        diff = norm(un-u(:,k))/norm(un);
        
        % UPDATE
        u(:,k) = ~channel_mask(:,k).*un + channel_mask(:,k).*u(:,k);
        
        % WRITE ON log FILE
        fprintf(fileID,'Channel %d, iter = %d, normalised difference of u: %2.4e\n',k,iter,diff);
        
        % TEST EXIT CONDITION
        if diff<tol
            break
        end
        
    end
    
end

fclose(fileID);

%% ---------------------------------------------------- GET THE 2D SOLUTION
u_end = reshape(u,iminfo.Height,iminfo.Width,colors);

%% ---------------------------------------------------- WRITE IMAGE OUTPUTS
imwrite(u_start,'masked_transport.png')
imwrite(u_end,'output_transport.png')

return

%% ------------------------------ AUXILIARY FUNCTION: ANISOTROPIC DIFFUSION
function u = anisodiff(u,dt,eps,geps,N,matrices)

% %% CHOICE 1: 
% M. Bertalmio, "Processing of flat and non-flat image information on 
% arbitrary manifolds using Partial Differential Equations", PhD Thesis, 2001.
for i=1:N
    ux = matrices.Dic*u;
    uy = matrices.Djc*u;
    uxx = matrices.M2i*u;
    uyy = matrices.M2j*u;
    uxy = matrices.M2ij*u;
    squared_normgrad = ux.^2 + uy.^2 + eps;
    u = u + dt*geps.*(uyy.*ux.^2 + uxx.* uy.^2 - 2*ux.*uy.*uxy)./squared_normgrad;
end

% %% CHOICE 2: 
% P. Perona, J. Malik, "Scale-space and edge detection using anisotropic 
% diffusion", PAMI 12(7), pp. 629-639, 1990.
% gamma = 1;
% for i=1:N
%     % image gradients in NSEW direction
%     uN=[u(1,:); u(1:m-1,:)]-u;
%     uS=[u(2:m,:); u(m,:)]-u;
%     uE=[u(:,2:n) u(:,n)]-u;
%     uW=[u(:,1) u(:,1:n-1)]-u;
% 
%     cN=1./(1+(abs(uN)/gamma).^2);
%     cS=1./(1+(abs(uS)/gamma).^2);
%     cE=1./(1+(abs(uE)/gamma).^2);
%     cW=1./(1+(abs(uW)/gamma).^2);
% 
%     u=u+dt*geps.*(cN.*uN + cS.*uS + cE.*uE + cW.*uW);

return