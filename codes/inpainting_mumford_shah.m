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


function [u_end,u_start]= inpainting_mumford_shah(imagefilename,maskfilename,maxiter,tol,param)
% Inpainting with the Mumford-Shah image model and Ambrosio-Tortorelli.
% For a given grey value image ustart with image domain \Omega and
% inpainting  domain (damaged part) D we want to reconstruct an image u from f
% by solving
% %%%%%%%%%%%%%% MINIMISATION PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = argmin_u  (\frac{\alpha}{2} \int_\Omega \chi^2 |\nabla u|^2 ~ dx %
%               + \beta \int_\Omega \left(\epsilon |\nabla \chi|^2     %
%               + \frac{(1-\chi)^2}{4\epsilon} \right)~ dx             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The above minimisation problem is solved iteratively via alternating
% solutions of the Euler-Lagrange equations for u and \chi.

%% ----------------- CREATE A log FILE WHERE TO STORE RESULTS IN txt FORMAT
logfilename = 'log_mumford_shah.log';
if exist(logfilename,'file')
    delete(logfilename);
end
fileID = fopen(logfilename,'w');

%% ------------------------------------ IMPORT THE CLEAN INPUT AND THE MASK
iminfo = imfinfo(imagefilename);
input  = im2double(imread(imagefilename));
% check if grayscale/truecolor dimension of image grey/colour
colors = size(input,3);

mask = im2double(imread(maskfilename));
mask = double(mat2gray(mask)==0); % characteristic function for the intact part of the image
if size(mask,3)==1 && colors>1
    mask = repmat(mask,[1,1,colors]);
end

%% ---------------------------------------------- GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;
N = iminfo.Height*iminfo.Width; % number of pixels

%% --------------------------------------------------------------- GRADIENT
% FORWARD AND BACKWARD
d1i_forward  = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[0,1],iminfo.Height,iminfo.Height)/h1;
d1j_forward  = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[0,1],iminfo.Width,iminfo.Width)/h2;
d1i_backward = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[-1,0],iminfo.Height,iminfo.Height)/h1;
d1j_backward = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[-1,0],iminfo.Width,iminfo.Width)/h2;
% PERIODIC BOUNDARY CONDITIONS
d1i_forward(end,:) = 0;
d1j_forward(end,:) = 0;
d1i_forward(end,[1 end]) = [1 -1]/h1;
d1j_forward(end,[1 end]) = [1 -1]/h2;
d1i_backward(1,:) = 0;
d1j_backward(1,:) = 0;
d1i_backward(1,[1 end]) = [1 -1]/h1;
d1j_backward(1,[1 end]) = [1 -1]/h2;

matrices.Dif  = kron(speye(iminfo.Width),d1i_forward);
matrices.Djf  = kron(d1j_forward,speye(iminfo.Height));
matrices.Dib  = kron(speye(iminfo.Width),d1i_backward);
matrices.Djb  = kron(d1j_backward,speye(iminfo.Height));

% CENTRAL
matrices.Dic = (matrices.Dif+matrices.Dib)/2;
matrices.Djc = (matrices.Djf+matrices.Djb)/2;

%% -------------------------------------------------------------- LAPLACIAN
d2i = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Height))/h1^2;
d2j = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Width))/h2^2;
% PERIODIC BOUNDARY CONDITIONS
d2i(1,end) = 1/h1^2;
d2i(end,1) = 1/h1^2;
d2j(end,1) = 1/h2^2;
d2j(1,end) = 1/h2^2;
% 2D domain LAPLACIAN
matrices.L = (kron(speye(iminfo.Width),d2i)+kron(d2j,speye(iminfo.Height)));

%% ------------------------------------------------------------ FREE MEMORY
clear d1i_forward dji_forward d1i_backward dji_backward
clear d2i d2j

%% -------------------------------------------------------------- ALGORITHM
% (1) INITIALIZATION: u^0 = 0, z^0 = 0, solve for k=1,2,...
u_start      = (mask).*input + ~mask;
channel_mask = reshape(mask,iminfo.Height*iminfo.Width,colors);
u            = reshape(u_start,iminfo.Height*iminfo.Width,colors);
chi          = zeros(N,colors);

rhsB = param.lambda*(channel_mask/param.alpha.*u);
rhsA = ones(N,colors); 

% FOR EACH COLOR CHANNEL
for k=1:colors
    
    % ITERATION
    for iter = 1:maxiter
        
        % SOLVE EULER-LAGRANGE EQUATION FOR u: 
        % i.e. B(\chi^k,u^k)= \chi_{\Omega\setminus D}/\alpha \cdot ustart.
        % B is a linear operator acting on u and reads
        % B(\chi,.) = -\mathrm{div}(\chi_\epsilon^2\nabla)
        %             + \chi_{\Omega\setminus D}/\alpha.
        % Solved via inversion of the linear operators.
        B      = matrixB(param,N,channel_mask(:,k),matrices,chi(:,k));
        unew   = B\rhsB(:,k);
        diff_u = norm(unew-u(:,k))/norm(unew);
        u(:,k) = unew;
        
        % SOLVE EULER-LAGRANGE EQUATION FOR \chi
        % i.e A(u^{k-1},\chi^k) = 1.
        % A is a linear operator acting on \chi and reads
        % A(u,.) = 1+\frac{2\epsilon\alpha}{\beta} |\nabla u|^2
        %           - 4\epsilon^2\Delta.
        % Solved via inversion of the linear operators.
        A        = matrixA(param,N,matrices,u);
        chinew   = A\rhsA(:,k);
        diff_chi = norm(chinew-chi(:,k))/norm(chinew);
        chi(:,k) = chinew;
        
        % WRITE ON log FILE
        fprintf(fileID,'Channel %d, normalised difference of u: %2.4e, normalised difference of chi: %2.4e\n',k,diff_u,diff_chi);
        
        % TEST EXIT CONDITION
        if diff_u<tol
            break
        end
    end    
end

fclose(fileID);

%% ---------------------------------------------------- GET THE 2D SOLUTION
u_end   = reshape(u,iminfo.Height,iminfo.Width,colors);
chi_end = mat2gray(reshape(chi,iminfo.Height,iminfo.Width,colors));

%% ---------------------------------------------------- WRITE IMAGE OUTPUTS
imwrite(u_start,'masked_mumford_shah_199.png')
imwrite(u_end,'output_mumford_shah_199.png')
% imwrite(chi_end,'levels_mumford_shah.png')

return

%% ---------------------------------------------------- AUXILIARY FUNCTIONS
function A = matrixA(param,N,matrices,u)
% Definition of (\nabla u)^2:
nablau2 = (matrices.Dic*u).^2 + (matrices.Djc*u).^2;

A = speye(N,N)...
    +2*param.epsilon*param.alpha/param.beta * spdiags(nablau2,0,N,N)...
    - 4*param.epsilon^2*matrices.L;
return

function B = matrixB(param,N,mask,matrices,chi)
% Definition of the nonlinear diffusion weighted by \chi^2:
a  = chi.^2+param.epsilon; % coefficient of nonlinear diffusion
ax = matrices.Dic*a;
ay = matrices.Djc*a;
A  = spdiags(a,0,N,N);
Ax = spdiags(ax,0,N,N);
Ay = spdiags(ay,0,N,N);

NonlinearDelta = A*matrices.L + Ax*matrices.Dic + Ay*matrices.Djc;

B = - NonlinearDelta + param.lambda*spdiags(mask,0,N,N)/param.alpha;
return


