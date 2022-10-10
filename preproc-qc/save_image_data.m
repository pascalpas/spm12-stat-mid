function save_nifti_file(input_matrix,v_b,output_filename,description)

%-----------------------------------------------------------------------------
% save_nifti_file.m
%
% Save-nifti_file will save any matrix in matlab workspace as an nifti file
%
% TO SAVE THE MATRIX AS ANALYZE, SET .NII TO .IMG IN FILENAME
%
% USAGE:
% 
% type : save_nifti_file(input_matrix,v_b,output_filename,description);
% 
% You are prompted to enter the:
% 1. name of the matrix you want to save
% 2. description of the data/matrix
% 3. target filename and location
%
% IMPORTANT: THE 'V_B' STRUCTURE SHOULD BE IN THE MATLAB WORKSPACE. THIS V_B 
% STRUCTURE SHOULD MATCH THE DATA MATRIX YOU WANT TO SAVE
%
% THIS 'V_B' STRUCTURE IS AUTOMATICALLY PRESENT IF YOU HAVE USED
% 'LOAD_nifti_FILES' OR 'GET_nifti_FILE'
%
% (c) Matthijs Vink, 10-2007
% version 1.0
%% -----------------------------------------------------------------------------
 
 warning off
  
  
  
if nargin<2 % if no input_matrix and/or v_b argument is given
disp('<save_nifti_file> Error! Specify which matrix to save!')
disp('<save_nifti_file> Error! Specify v_b structure!')
disp('')
disp('<save_nifti_file> Type : save_nifti_file(input_matrix,v_b)')
disp('<save_nifti_file> or help save_nifti_file ')

elseif nargin<4 % if not all input arguments given

  % if input_nifti_files is unspecified
  [description] = input('<save_nifti_file> Enter a description of the data : ','s');
   % determine the file location and name
  [filename2, pathname2] = uiputfile('*.nii', 'Save the contrast');
  output_filename = [pathname2,filename2];

 
disp('<save_nifti_file> Assigning new filename and file description')

v_b.fname=[output_filename];
v_b.descrip =[description];

disp('<save_nifti_file> Assigning matrix dimensions')

v_b.dim(1,1:3)=[size(input_matrix,1) size(input_matrix,2) size(input_matrix,3)];

disp('<save_nifti_file> Writing nifti file')

spm_write_vol(v_b,input_matrix);

disp('<save_nifti_file> Done!')

elseif nargin>3 
 
disp('<save_nifti_file> Assigning new filename and file description')

v_b.fname=[output_filename];
v_b.descrip =[description];

disp('<save_nifti_file> Assigning matrix dimensions')

v_b.dim(1,1:3)=[size(input_matrix,1) size(input_matrix,2) size(input_matrix,3)];

disp('<save_nifti_file> Writing nifti file')

spm_write_vol(v_b,input_matrix);

disp('<save_nifti_file> Done!')

end

