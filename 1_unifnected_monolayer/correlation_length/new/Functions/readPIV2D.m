function [x,y,u,v,s2n]=readPIV2D(filename,outl,s2nl,umax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READPIV2D - reads 2D PIV vectors
%           - filters outlayers
%           - checks signal-to-noise ratio  
%
%    Inputs:
%
%    filename - name of the file to be read
%    outl     - flag, if 1 outlayers are removed
%                     if 0 outlayers are not removed
%    s2nl     - flag, if >0 points with signal-to-noise < s2nl are removed
%                     if =0 signal-to-noise is not checked
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %
  % loads data
  %
  vec = load(filename);   

  %
  % determines sizes of the 3D vector array
  % 
  n1 = length(unique(vec(:,1)));
  n2 = length(unique(vec(:,2)));

  %
  % reshapes data into 3D format and allocates
  % into different variables 
  %
  vec=reshape(vec,n1,n2,5);
  
  x    = vec(:,:,1);
  y    = vec(:,:,2);
  u    = vec(:,:,3);
  v    = vec(:,:,4);
  s2n  = vec(:,:,5);

  %
  % s2n correction
  %
  if s2nl>0

  % Remove bad s2n points

     u(s2n<s2nl)=NaN;
     v(s2n<s2nl)=NaN;

     uu=sqrt(u.^2+v.^2);
     u(uu>umax)=NaN;
     v(uu>umax)=NaN;
     %for iz=1:n3
     %   u(:,:,iz) = inpaint_nans(u(:,:,iz),3);
     %   v(:,:,iz) = inpaint_nans(v(:,:,iz),3);
     %   w(:,:,iz) = inpaint_nans(w(:,:,iz),3);
     %end
    
  end 
  %
  % Adaptive Local Median filtering
  %
  % smof for 1000*1000 image is 5.0000e-05
  % Z = SMOOTHN(Y,S) smoothes the array Y using the smoothing parameter S.
  %   S must be a real positive scalar. The larger S is, the smoother the
  %   output will be.
  if outl>0
     %u = medfilt2(u,[3 3]);
     %v = medfilt2(v,[3 3]); 
     %%
     smof = 0.05/sqrt(size(u,1)*size(u,2));
  % how to choose the smof.
     u = smoothn(u,smof,'robust','MaxIter',100,'TolZ',1e-3);
     v = smoothn(v,smof,'robust','MaxIter',100,'TolZ',1e-3);
  end
   
