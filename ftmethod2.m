function [pangle,I,In,thetarad,F, Ipl]=ftmethod2(f)
% [I, IN, THETARAD, F] = FTMETHOD applies the fourier transfrom method to a
% grayscale image f. If no image is supplied in the argument then the user 
% is prompted to select a file from a gui. The function returns I, the
% intensity, the area of I, THETARAD, the associated angles, and F, the orientation
% tensor. It is advisable to crop your image so that it has the same number
% of pixels in the vertical and horizontal directions, i.e. a square. This
% has to do with the fact that the frequencies are dependent on the image
% size.
%
% ADDTIONAL FUNCTIONS - The following functions are custom 
%
%   gscale.m - performs contrast stretching
%   winfilt.m - creates a windowing filter, in this case a Welch window
%   hpfilter.m - creates a highpass filter
%   lpfilter.m - creates a lowpass filter
%   pfft.m  - converts the magnitude spectrum from a caretsion
%             representation to a polar one expressed in matrix form.
%   integrate.m - applies Simpson's 1/3 rule to estimate the integral
%
%
% DETAILS - The imported image first undergoes constrast stretching to set
% the lowest value in the image to 0 and the highest value to 255. The image
% is then filtered with a Welch window to eliminate spectral leakage in the
% FFT. The windowed image is then converted to the frequency space with
% FFT2 and filtered with a low and high pass filter. The filtered fft is
% then transformed to polar coordinates and binned according to the spacing
% of r and theta set in pfft.m into the matrix h. A ones matrix the same 
% size as f is filtered with the same high and low pass filters and 
% converted to a polar representation. The resulting matrix, hn, contains 
% the number of points that contributed to the value in the bin for a given
% r and theta. Dividing h by hn gives an average intensity. Finally, the
% intensity values along r for a given theta are summed and divided by the
% number of r values for that theta. These values are offset 90 degrees in
% the frequency space from the orientation in the spatial domain. So theta
% is rotated 90 degrees, thus giving I and thetarad, which is the
% orientation distribution for the image. The orientation matrix is then
% determined from weighted integrals of I(thetarad).
%
% Written by Dr. Edward A Sander, Department of Biomedical Engineering
% University of Minnesota, Minneapolis, MN.
%
% For more details, see:
%
% Sander, E.A. and Barocas, V.H. 2008, "Comparison of 2D Fiber Network
% Orientation Measurement Methods," Journal of Biomedical Materials
% Research Part A, 88, pp. 322-331.


%%%%%%%% Welch Window is ON !!!!!!!!!!!!!

if nargin == 0
    [filename, pathname] = uigetfile('*.*','Choose Image to Analyze');
    f = imread(fullfile(pathname,filename));
end

f = f(:,:,1);


%%% Contrast stretch the image

f=gscale(f);
f = double(f);
[M,N]=size(f);

%%%%%%%%%%%% NoT in Use Now %%%%%%%%%%%%%%
% Make the image even in size
% if mod(M,2)==1
%     f=f(2:end,:);
% end
% 
% if mod(N,2)==1
%     f=f(:,2:end);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Use Welch Window on Image
% Note: Comment out one of the FF lines immediately below depending on
% whether you want the window on or not.
W = winfilt(f,'welch');
FF = fft2(W);
% FF = fft2(f);
[M,N] = size(f);

% Set up the cutoff frequencies for the filters
if M > N
    lowcutoff = round(0.0156*M); % This factor 0.0156 is somewhat arbitrary
    highcutoff = round(M/4);
else
    lowcutoff = round(0.0156*N);
    highcutoff = round(N/4);
end
% NOTE: The cutoff frequencies are dependent on the image size and the 
% object size. They can be determined as freq=ImageSize/(2*fiber diameter).
% For example, if you have an image that is 512x512 and you are interested 
% in fibers that are 10 pixels in diameter. Then freq = 512/(2*10)= 25.6.
% So make sure that your cutoff frequencies do not exclude this value.

Hh = hpfilter('ideal',M,N,lowcutoff);
Hl = lpfilter('ideal',M,N,highcutoff);

FF2 = FF.*Hh.*Hl;
K = ones(M,N);
K2 = K.*Hh.*Hl;

[h,theta,rho] = pfft(FF2);
[hn]= pfft(K2);
%figure, imshow(ifft2(FF2));
% If you want to see Fourier domain
%figure, imshow(log(1+abs(fftshift(FF2))),[])

% Find how many points contribute to each theta
[x,y] = find(hn);
val = ones(size(x));
[RR,TT]=size(hn);
NN = full(sparse(x,y,val,RR,TT));
% As you move out radially, more points become associated with the line.
% This has the effect of including more higher than lower frequency 
% values. To remedy this, we divide the cumulutive intensity value
% contained in h by the number of points that contributed to that value
% contained in hn. We then sum these values and divide by the number of r
% occupied r bins (NN), which may vary with theta.

I = sum(h./(hn+eps))./sum(NN+eps); %Note add the eps to avoid dividing by 
% zero. Where the denominator is zero the numerator should also be zero,
% thus avoiding the introduction of very high terms.
% I = sum(h)./sum(hn+eps); %This divides the sum of the amplitudes by the
% number of points on r. 

%%%% Orientation Matrix %%%%%

[r,c]=find(theta==0);% Need to rotate 90 degrees because the FT is 90 off but 
I=[I(1,c:end) I(1,2:c)];%keep the integration domain between -pi/2 and pi/2.
thetarad = theta*pi/180;%
CosTop = integrate(I.*cos(thetarad).^2,thetarad);
SinTop = integrate(I.*sin(thetarad).^2,thetarad);
CSTop = integrate(I.*sin(thetarad).*cos(thetarad),thetarad);
In = integrate(I,thetarad);

F = [CosTop CSTop;CSTop SinTop]./In;
[FV,FD]=eig(F);

if FD(1) > FD(4)
    pangle = atan(FV(2)/FV(1))*180/pi;
    if sign(FV(2))+sign(FV(1)) == 0
        if pangle > 0
            pangle = pangle+90;
        else
            pangle = pangle+180;
        end
    end
else
    pangle = atan(FV(4)/FV(3))*180/pi;
    if sign(FV(4))+sign(FV(3)) == 0
        if pangle >0
            pangle = pangle + 90;
        else
            pangle = pangle + 180;
        end
    end
end
%pangle;
%%%% Plot vectors on picture
%{
figure(1), imshow(f,[])
hold on
p=[N/2+1;N/2+1];
q=[M/2+1; M/2+1];
u=[FD(1)*FV(1);FD(4)*FV(3)];
v=[FD(1)*FV(2);FD(4)*FV(4)];
quiver(p,q,u,-v,0.5*M,'LineWidth',3,'Color','r'); % Note have to flip the magnitude of y because of 
% the way the image is indexed - increasing y points down in the image.
hold off
%}

%%% Generate polar plot of distribution

Ipl = [I./(In) I(2:end)./(In)];

tpl = [thetarad (thetarad(2:end)+pi)];
%figure, polar(tpl,Ipl) 

%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%
function H = lpfilter(type, M, N, D0, n)
%LPFILTER Computes frequency domain lowpass filters.
%   H = LPFILTER(TYPE, M, N, D0, n) creates the transfer function of
%   a lowpass filter, H, of the specified TYPE and size (M-by-N). To
%   view the filter as an image or mesh plot, it should be centered
%   using H = fftshift(H). 
%
%   Valid values for TYPE, D0, and n are:
%
%   'ideal'    Ideal lowpass filter with cutoff frequency D0. n need
%              not be supplied.  D0 must be positive.
%
%   'btw'      Butterworth lowpass filter of order n, and cutoff
%              D0.  The default value for n is 1.0.  D0 must be
%              positive.
%
%   'gaussian' Gaussian lowpass filter with cutoff (standard
%              deviation) D0.  n need not be supplied.  D0 must be
%              positive. 

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.7 $  $Date: 2003/10/13 00:46:25 $

% Use function dftuv to set up the meshgrid arrays needed for
% computing the required distances. 
[U, V] = dftuv(M, N);

% Compute the distances D(U, V).
D = sqrt(U.^2 + V.^2);

% Begin filter computations.
switch type
case 'ideal'
   H = double(D <= D0);
case 'btw'
   if nargin == 4
      n = 1;	
   end
   H = 1./(1 + (D./D0).^(2*n));
case 'gaussian'
   H = exp(-(D.^2)./(2*(D0^2)));
otherwise
   error('Unknown filter type.')
end


function H = hpfilter(type, M, N, D0, n)
%HPFILTER Computes frequency domain highpass filters.
%   H = HPFILTER(TYPE, M, N, D0, n) creates the transfer function of a
%   highpass filter, H, of the specified TYPE and size (M-by-N). 
%   Valid values for TYPE, D0, and n are:
%
%   'ideal'     Ideal highpass filter with cutoff frequency D0. n need not
%               be supplied. D0 must be positive.
%
%   'btw'       Butterworth highpass filter with cutoff (standard 
%               deviation) D0. The default value for n is 1.0. D0 must be 
%               positive.
%
%   'gaussian'  Gaussian highpass filter with cutoff (standard deviation)
%               D0. n need not be supplied. D0 must be positive.
%
%

% The transfer function Hhp of a highpass filter is 1-Hlp, where Hlp is the
% transfer function of the corresponding lowpass filter. Thus, we can use
% function lpfilter to generate highpass filters.

if nargin ==4
    n = 1; % Default value of n.
end

% Generate highpass filter.
Hlp = lpfilter(type,M,N,D0,n);
H = 1- Hlp;

function g = gscale(f, varargin)
%GSCALE Scales the intensity of the input image
%   G = GSCALE(F,'full8') scales the intensities of F to the full 8-bit
%   intensity range [0, 255]. This is the default if there is only one
%   input argument.
%
%   G = GSCALE(F,'minmax', LOW, HIGH) scales the intensities of F to the
%   range [LOW, HIGH]. These values must be provided, and they must be in
%   the range [0,1], independently of the class of the input. GSCALE
%   performs any necessary scaling. If the input is of class double, and
%   its values are not in the range [0, 1], then GSCALE scales it to this
%   range before processsing.

if length(varargin) == 0
    method = 'full8';
else
    method = varargin{1};
end

if strcmp(class(f), 'double') & (max(f(:)) > 1 | min(f(:)) < 0)
    f=mat2gray(f);
end

% Perform the specialized scaling.
switch method
    case 'full8'
        g = im2uint8(mat2gray(double(f)));
    case 'full16'
        g = im2uint16(mat2gray(double(f)));
    case 'minmax'
        low = varargin{2}; high = varargin{3};
        if low > 1 | low < 0 | high > 1 | high < 0
            error('Parameters low and high must be in the range [0, 1].')
        end
        if strcmp(class(f),'double')
            low_in = min(f(:));
            high_in = max(f(:));
        elseif strcmp(class(f),'uint8')
            low_in = double(min(f(:)))./255;
            high_in = double(max(f(:)))./255;
        elseif strcmp(class(f),'uint16')
            low_in = double(min(f(:)))./65535;
            high_in = double(max(f(:)))./65535;
        end
        % imadjust automatically matches the class of the input. 
        g = imadjust(f, [low_in high_in], [low high]);
    otherwise
        error('Unknown method.')
end

function [I]=integrate(f,x,type)
%INTEGRATE takes the row vectors F and X and applies one of the following numerical
%integration schemes:
%
% See Chapra and Canale Numerical Methods for Engineers 2002.
%   'SIMP3'  -   Simpson's 1/3 Rule. Default Method.
%   'TRAP'  -   Trapezoidal Rule. 

f=f(:)'; % Make sure that f and x are row vectors
x=x(:)';

if size(f) ~= size(x)
    error('vectors f and x are not the same size')
end

% Default integration scheme.
if nargin < 3
    type = 'SIMP3';
end


% General Variables

ind = size(f,2);
h = (x(end) - x(1))/(ind - 1);

switch type
    case 'SIMP3'
        fodd = f(3:2:end-1);
        feven = f(2:2:end-1);
        I = h/3*(f(1) + 4*sum(feven) + 2*sum(fodd) + f(end));
        
    case 'TRAP'
        I = 0.5*h*(f(1) + 2*(sum(f(2:end-1)) + f(end)));
    otherwise
        error('Unknown scheme')
end

function [h,theta,rho]= pfft(f)
%PFFT - Transforms the fourier magnitude spectrum to polar coordinates.
% Code adapted from hough.m in Digital Image Processing Using Matlab -
% Gonazalez et al.


drho = 1;
dtheta = 1;
%%%%%%%%%%%%%%%%%%
%     1   N/2   N
%  1  -----------
%     |  I | II |   I is fbot and contains the dc component at 1,1
% M/2 -----------   IV is ftop. Reshuffling the matrix is the same as using
%     | IV | III|   fftshift.
%  M  -----------
%%%%%%%%%%%%%%%%%%%

%%%%%%%  Reminder on how the angle is calculated here.
%+90       Q is theta
%|   /
%|  /
%| /
%|/ Q)
%-----------0
%|\ Q)
%| \
%|  \
%|   \ 
%-90
% Break up the Fourier space into the I and IV quadrants
f=abs(f).^2; % Get power spectrum
[M,N] = size(f);
ftop = f(ceil(M/2)+1:end,1:ceil(N/2));
fbot = f(1:ceil(M/2),1:ceil(N/2)); % Note this is same as 0 to M/2-1 

% Set up the matrix of rho and theta values
theta = linspace(-90,0,ceil(90/dtheta) +1);
theta = [theta -fliplr(theta(1:end-1))];
ntheta = length(theta);
D = sqrt((M/2-1)^2+(N/2-1)^2); % I'm wondering if I should subtract 0.5 instead of -1.
q = ceil(D/drho);
nrho = q+1;
rho = linspace(0, q*drho, nrho);

% Let's just work on fbot
[yb, xb, valb] = find(fbot); % row, col
xb = xb - 1; yb = yb - 1;
rho_matrix_b = sqrt(xb.^2+yb.^2);
theta_matrix_b = - atan(yb./(xb+eps))*180/pi;

% Now ftop
fl_ftop =flipud(ftop);
[yt, xt, valt] = find(fl_ftop);
xt = xt-1;
rho_matrix_t = sqrt(xt.^2+yt.^2);
theta_matrix_t = atan(yt./(xt+eps))*180/pi;

% Combine

rho_matrix = [rho_matrix_t; rho_matrix_b];
theta_matrix = [theta_matrix_t; theta_matrix_b];
val = [valt; valb];

% Initialize
h = zeros(nrho, length(theta));

rho_slope=(nrho-1)/(rho(end)-rho(1));
rho_bin_index = round(rho_slope*(rho_matrix-rho(1))+1);
theta_slope = (ntheta-1)/(theta(end)-theta(1));
theta_bin_index = round(theta_slope*(theta_matrix-theta(1))+1);

h = h + full(sparse(rho_bin_index(:),theta_bin_index(:),...
    val(:),nrho,ntheta));

%imshow(log(1+h),[],'XData',theta,'YData',rho,'InitialMag','fit')
%axis on, axis normal

function [h]=winfilt(f,type)
% WINFILT produces a windowing function for the image f. 
%
%
% Valid window funcitons for TYPE are:
%
% 'hann'    Hanning Window - this is default TYPE
%
% 'welch'   Welch Window

if nargin == 1
    type = 'hann'
end

[M,N] = size(f);
Tr = M-1;
Tc = N-1;

r = 0:Tr;
c = 0:Tc;


switch type
    case 'hann'
        mr = 0.5 - 0.5*cos(2*pi*r/Tr);
        nc = 0.5 - 0.5*cos(2*pi*c/Tc);
    case 'welch'
        mr = 1 - ((r-0.5*Tr)./(0.5*Tr)).^2;
        nc = 1 - ((c-0.5*Tc)./(0.5*Tc)).^2;
    otherwise
        error('Unknown filter type')
end

        
R = repmat(mr',1,N);
C = repmat(nc,M,1);

H = R.*C;
h = H.*f;

function [U, V] = dftuv(M,N)
%DFTUV Computes meshgrid frequency matrices.
%   [U, V] = DFTUV(M, N) computes meshgrid frequency matrices U and V. U
%   and V are useful for computing frequency-domain filter functions that
%   can be used with DFTFILT. U and V are both M-by-N.

% Set up range of variables.
u = 0:(M-1);
v = 0:(N-1); 

% Compute the indices for use in meshgrid. 
idx = find(u > M/2);
u(idx) = u(idx) - M;
idy = find(v > N/2);
v(idy) = v(idy) - N;

% Compute the meshgrid arrays.
[V, U] = meshgrid(v, u);
