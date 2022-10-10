function [images,v_b,xyz_b,y_b] = get_image_data2(filename,commandline)
% -----------------------------------------------------------------------------
% get_image_data2.m
%
% get_image_data2 is a function and will load an NIFTI file and place the data into matlab
% workspace
%
% USAGE:
% 
% type : [images,v_b,xyz_b,y_b] = get_image_data2(filename,message);
%
% You are prompted to select the 3d nifti file you want to load
%
% 
% example: [images,v_b,xyz_b,y_b] = get_image_data2([],message);
%          which will prompt a gui with a specified message
% 
% example: [images,v_b,xyz_b,y_b] = get_image_data2(files,message);
%          which will load the specified files, and ignores message
% 
% example: [images,v_b,xyz_b,y_b] = get_image_data2(files);
%          which will load the specified files
% 
%
% The following variables will appear in the matlab workspace:
%
% images : complete path and filename of selected NIFTI file
% v_b    : SPM5 structure with all relevant image information (eg voxel size)
% xyz_b  : SPM5 structure with all voxel coordinates 
% y_b    : the actual data matrix
%
%
% (c) Matthijs Vink, 01-2009
% version 1.0
% -----------------------------------------------------------------------------

if nargin < 1
[images,sts] = spm_select(inf,'image','Select Nifti file');
elseif nargin < 2 
    images = filename;
elseif nargin == 2
     if isempty(filename)==1
     [images,sts] = spm_select(inf,'image',[commandline]);
     else
      images = filename;   
     end    
end

v_b = spm_vol([images]);
[y_b,xyz_b] = spm_read_vols(v_b);
    

