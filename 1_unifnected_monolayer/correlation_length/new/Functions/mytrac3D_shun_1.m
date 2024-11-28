function [txx,tyy,tz,ut00,vt00,wt0,X0,Y0,Lx,Ly,elaser_tan, elaser_tann,elaser_norm, elaser_x, elaser_y]...
    = mytrac3D_shun_1(E,h,h0,nx0,ny0,fil,fcalx,Nyq)
%function for calculating traction stress on the surface with the boundary
%condition w=0; Applied to the experimental setup with agar on top.
%This is a psuedo 3D TFM
%Can Use modified DIC image to confine the deformation around the cell
%Don't need to -xr(1)/2 to put vectors in the center because window didn't
%crop






%E-Young's modulus
%h-thickness of gel
%h0-distance from the bottom to the beads layer, generally
%h0=h|h-0.5|h0-1,here it is the depth, not layers.
%nx0 256/128
%ny0 256/128
%fil-file name for the PIV results
%fcalx-scale in xy
%fcalz-scale in z
%xcg-x component of the center
%ycg
%II-length of the box
%wind-flag for filtering
%Nyq-2 for 50% overlapping, 4 for 75% overlapping
%mic
%bb-DIC image, use to get the contour

%tx-traction force in x
%ut0-velocity vetor in x before filtering
%X0-grid generated
%Lx-box size
%elaser_tan-tangential elastic energy(non-dimensional)



% se1  = strel('disk',     2);
% se2  = strel('disk',     10);

%[x,y,z,u,v,w,s2n]=readpiv3d(fil,0,0,5);
% [x,y,z,u,v,w,s2n]=readPIV3D(fil,1,1,2);
%%change these values depedning on how large the displacements are expected to be and what can be considered as outlier 
  [x,y,u,v,s2n]=readPIV2D(fil,1,1,5); 
 %%%% ADAPT HERE (Maybe outliers, signal to noise, remove displacements larger than 5 pixels.)
 %%%% If cells move more, increase last value, because all displaements higher than that are cut out
  
% 
%   im = bb>0;
%   im = double(im);
%   
  
  u(isnan(u))=0;
  u=u-median(u(:));
  v(isnan(v))=0;
  v=v-median(v(:));
%   
%   u = u.*im(16:8:984,16:8:984)';
%   v = v.*im(16:8:984,16:8:984)';
        
        
        
        
  
  n1= size(x,1);
  n2= size(y,1);
  
  % median
  u(isnan(u))=0;
  u=u-median(u(:));
  v(isnan(v))=0;
  v=v-median(v(:));
  
  w = zeros(size(u));

% %%%%%%%%%%%%

xr = sort(unique(x))*fcalx;
yr = sort(unique(y))*fcalx;


ix0 =  1; 
ixe = n1;
jy0 =  1; 
jye = n2; 

vecx=ix0:ixe;
vecy=jy0:jye;

%theory behind it
KXM = pi/(xr(2)-xr(1))*2/Nyq;
KYM = pi/(yr(2)-yr(1))*2/Nyq;

% %new xy has real length, x has the values 16,32,64... -1/2 to put in center
% x  = xr(vecx) - xr(1)/2;
% y  = yr(vecy) - yr(1)/2;

x = xr;
y = yr;



% u = u(vecx,vecy)*fcalx;
% v = v(vecx,vecy)*fcalx;

u = u.*fcalx;
v = v.*fcalx;


 
 Lx0 = (x(end) - x(1)); %bego
 Ly0 = (y(end) - y(1)); %bego

xgap = Lx0/(nx0-1);
ygap = Ly0/(ny0-1);

%X0 Y0 are in um
[Y0,X0] = meshgrid(y(1):ygap:y(end),x(1):xgap:x(end));
[y0,x0] = meshgrid(y,x);



ut = interp2(y0,x0,u,Y0,X0);
vt = interp2(y0,x0,v,Y0,X0);
wt = interp2(y0,x0,w,Y0,X0);

ut0=ut;
vt0=vt;
wt0=wt;

%       fudges for extrap

ut(ut~=ut)=0;
vt(vt~=vt)=0;
wt(wt~=wt)=0;


% X Y is the left shift for XO and Y0, X and XO has the same size.
[Y,X] = meshgrid(0:ygap:Ly0,0:xgap:Lx0);


% filtering
 %alpha = 0.2 ;%%%if you change alpha of tukey window it will decay differently to zero
 alpha = 0.1;
% 
HWx = ones(length(X')) ;

HWy=HWx;
% 
 HWx(X<=alpha*(Lx0)/2) = (1 + cos(pi*(2*X(X<=alpha*(Lx0)/2)'/(alpha*(Lx0))-1)))/2 ;
 HWx((X>alpha*(Lx0)/2).*(X<=(Lx0)*(1-alpha/2))>0) = 1 ;
 HWx(X>(Lx0)*(1-alpha/2)) = (1 + cos(pi*(2*X(X>(Lx0)*(1-alpha/2))'/(alpha*(Lx0))-2/alpha+1)))/2 ;
% 
 HWy(Y<=alpha*(Ly0)/2) = (1 + cos(pi*(2*Y(Y<=alpha*(Ly0)/2)'/(alpha*(Ly0))-1)))/2 ;
 HWy((Y>alpha*(Ly0)/2).*(Y<=(Ly0)*(1-alpha/2))>0) = 1 ;
 HWy(Y>(Ly0)*(1-alpha/2)) = (1 + cos(pi*(2*Y(Y>(Ly0)*(1-alpha/2))'/(alpha*(Ly0))-2/alpha+1)))/2 ;
% 
HW=sqrt(HWx.*HWy);
% figure;
% imagesc(HW)
% pause()

    ut = ut.*HW;
    vt = vt.*HW;
    wt = wt.*HW;
    %
    nx = nx0;
    ny = ny0;
    Lx = Lx0;
    Ly = Ly0;



    Theta_x = 2*pi*X'/Lx0;
    Theta_y = 2*pi*Y'/Ly0;
    
    %	Transforms displacements to Fourier space
%
ut =     fft2(ut);  %the zero-frequency component is in the upper-left corner of the two-dimensional FFT
% ut = fftshift(ut);
%
vt =     fft2(vt);  %the zero-frequency component is in the upper-left corner of the two-dimensional FFT
% vt = fftshift(vt);
%
wt =     fft2(wt);  %the zero-frequency component is in the upper-left corner of the two-dimensional FFT


ut = fftshift(ut);  %the zero-frequency component is near the center of the matrix
vt = fftshift(vt);  %the zero-frequency component is near the center of the matrix
wt = fftshift(wt);  %the zero-frequency component is near the center of the matrix

kx = [-nx/2:nx/2-1]*2*pi/Lx;
ky = [-ny/2:ny/2-1]*2*pi/Ly;
k2max = max(kx)^2 + max(ky)^2; 

KXM = pi/(xr(2)-xr(1))*2/Nyq;
KYM = pi/(yr(2)-yr(1))*2/Nyq;


nfil = 2; %width of the Gaussian filter in multiples of (dx^2 + dy^2)^0.5    
%	Compute the tractions from the displacements

sig = 0.46; %Poisson ratio, perfectly elastic material would be 0.5, how much width is reduced when stretched

I = sqrt(-1);
tx = zeros(nx,ny);
ty = zeros(nx,ny);
tz = zeros(nx,ny);




for ix=1:nx
   if (abs(kx(ix))<=KXM)
      for iy=1:ny
         if (abs(ky(iy))<=KYM)
            k2 = kx(ix)^2+ky(iy).^2; 
            filter = exp(nfil*k2/k2max/4*pi^2);
%            [txz,tyz,tzz,uu,uv,uw] = forcesonlytzz(kx(ix),ky(iy),h,sig,squeeze(tz(ix,iy)),E);
%            [txz,tyz,tzz] = forcestzz(kx(ix),ky(iy),h,h,sig,ut(ix,iy),vt(ix,iy));
           % [txz,tyz,tzz] = forces3D(kx(ix),ky(iy),h,h,sig,ut(ix,iy),vt(ix,iy),mic*wt(ix,iy));
           [txz,tyz,tzz] = forces3D(kx(ix),ky(iy),h,h0,sig,ut(ix,iy),vt(ix,iy),wt(ix,iy));
            tx(ix,iy) = E*txz/filter;
            ty(ix,iy) = E*tyz/filter;
            tz(ix,iy) = E*tzz/filter;
         end
      end
   end
end

% 
% zero modes
%
tx(nx/2+1,ny/2+1)=ut(nx/2+1,ny/2+1)*E/h/(1+sig)/2;
ty(nx/2+1,ny/2+1)=vt(nx/2+1,ny/2+1)*E/h/(1+sig)/2;
tz(nx/2+1,ny/2+1)=wt(nx/2+1,ny/2+1)*E*(1-sig)/h/(1+sig)/(1-2*sig);


%	Brings the tractions to physical space

tx = fftshift(tx);
ty = fftshift(ty);
tz = fftshift(tz);

% elaser=sum(sum(tx*ut+ty*vt+tz*wt));
% elaser=real(elaser);

tx = real(ifft2(tx));
ty = real(ifft2(ty));
tz = real(ifft2(tz));

ux = real(ifft2(ut));
uy = real(ifft2(vt));  
uz = real(ifft2(wt));

Lx =Lx0;%%%%%%%%%we need to put this box size to find comparable results, because we are taking the same number of points and the stress matrices have the same size.
Ly =Ly0;
% elaser=sum(sum(tx*ux+ty*uy+tz*uz));

  

%for energy do we need a 0.5?
elaser_tot=sum(sum(tx.*ut0+ty.*vt0+tz.*wt0));
txx=tx;tyy=ty;ut00=ut0;vt00=vt0;
elaser_tann = sum(sum(txx.*ut00+tyy.*vt00)); 
elaser_tan = sum(sum(tx.*ut0+ty.*vt0));
elaser_x = sum(sum(tx.*ut0));
elaser_y = sum(sum(ty.*vt0));
elaser_norm = sum(sum(tz.*wt0));