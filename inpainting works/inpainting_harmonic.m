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

function [u_end,u_start]=inpainting_harmonic(imagefilename,maskfilename,lambda,tol,maxiter,dt)

%% ----------------- CREATE A log FILE WHERE TO STORE RESULTS IN txt FORMAT
logfilename = 'log_harmonic.log';
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
mask = double(mat2gray(mask)>0);
if size(mask,3)==1 && colors>1
    mask = repmat(mask,[1,1,colors]);
end

%% ---------------------------------------------------- INITIALIZATION OF u
u_start = (~mask).*input + mask;
u_end   = zeros(iminfo.Height,iminfo.Width,colors);

%% ---------------------------------------------- GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1; 

%% -------------------------------------------------------------- LAPLACIAN
d2i = toeplitz(sparse([1,1],[1,2],[-2,1]/h1^2,1,iminfo.Height));
d2j = toeplitz(sparse([1,1],[1,2],[-2,1]/h2^2,1,iminfo.Width));
% NEUMANN BOUNDARY CONDITIONS
d2i(1,[1 2])         = [-1 1]/h1;
d2i(end,[end-1 end]) = [1 -1]/h1;
d2j(1,[1 2])         = [-1 1]/h2;
d2j(end,[end-1 end]) = [1 -1]/h2;
% 2D domain LAPLACIAN
L = kron(speye(iminfo.Width),d2i)+kron(d2j,speye(iminfo.Height));

%% ------------------------------------------------------------ FREE MEMORY
clear d2i d2j

%% -------------------------------------------------------------- ALGORITHM
% INITIALIZATION
u            = reshape(u_start,iminfo.Height*iminfo.Width,colors);
f            = reshape(u_start,iminfo.Height*iminfo.Width,colors);
channel_mask = reshape(mask,iminfo.Height*iminfo.Width,colors);

% FOR EACH COLOR CHANNEL
for k=1:colors
    for iter = 1:maxiter
        % COMPUTE NEW SOLUTION
        unew = u(:,k) + dt*(L*u(:,k) + lambda*(~channel_mask(:,k)).*(f(:,k)-u(:,k)));
        
        % COMPUTE EXIT CONDITION
        diff = norm(unew-u(:,k))/norm(unew);
        
        % UPDATE
        u(:,k) = unew;
        
        % TEST EXIT CONDITION
        if diff<tol
            break
        end
    end
    
    % WRITE ON log FILE
    fprintf(fileID,'Channel %d: Iterations: %d, Normalised difference of u: %2.4e\n',k,iter,diff);
    u_end(:,:,k) = reshape(u(:,k),iminfo.Height,iminfo.Width);
end

fclose(fileID);

%%  ---------------------------------------------------- WRITE IMAGE OUTPUT
imwrite(u_start,'masked_harmonic.png')
imwrite(u_end,'output_harmonic.png')

return