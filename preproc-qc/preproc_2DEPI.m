function preproc_2DEPI(subjectDir,subjectName,anatomy,func,scanTR,mbFactor)
% Input:
%   subjectDir:     directoryname of subject;
%   subjectName:    subject name;
%   anatomy:        anatomical dir name
%   func:           functional dir name
%   scanTR:         TR of functionals (e.g. 1)
%   mbFactor:       multiband factor for slice timing  (e.g. 3)

% e.g.: preproc_2DEPI('nfs/home2/adima/ses-1002500B_share','Juan','anat','func',3,2)


warning('off');

anatfile = dir(fullfile(subjectDir,anatomy,['*.nii']));
funcfile = dir(fullfile(subjectDir,func,['*.nii']));
[~,scanCode,~] = fileparts(funcfile.name);
[~,anatCode,~] = fileparts(anatfile.name);

hp_cutoff = 128;

if length(funcfile) > 1 
    funcHdr = NaN;
    funcData = NaN;
    maxval = NaN;
else
    funcHdr = spm_vol( fullfile(subjectDir,func,funcfile.name) );
    funcData = spm_read_vols( funcHdr );
    maxval = max(max(max(max(funcData))));
    
    disp('>> >> Converting 4D to 3D nifti files')
    disp(' ');
    spm_file_split(fullfile(subjectDir,func,funcfile.name));
    eval(['!mv ',fullfile(subjectDir,func,funcfile.name),' ',fullfile(subjectDir,func,['_',funcfile.name])]);
    
    %remove first two dynamics
     toCopy = dir(fullfile(subjectDir,func,['*0003.nii']));
     scan3 = toCopy.name;
     scan1 = [scan3(1:end-5) '1.nii'];
     scan2 = [scan3(1:end-5) '2.nii'];
     eval(['!cp ',fullfile(subjectDir,func,scan3),' ',fullfile(subjectDir,func,scan1)])
     eval(['!cp ',fullfile(subjectDir,func,scan3),' ',fullfile(subjectDir,func,scan2)])
end

outputdir = fullfile(subjectDir,'qc_report');
if exist(outputdir,'dir') ~= 7
    mkdir (outputdir)
end

%% PREPROC %%
clear matlabbatch
localPath = which('preproc_2DEPI.m');
[localDir,~,~] = fileparts(localPath);
sample_jobman = fullfile(localDir,'sample_jobman_preproc_SPM12_2d_epi.mat');
load(sample_jobman);

% Change files
% Input for slice time correction
funcpath = fullfile(subjectDir,func);
anatpath = fullfile(subjectDir,anatomy);
eval(['[input_slicetime,dirs]=spm_select(''FPList'',funcpath,''^',scanCode,'.*\.nii$'');']);

[images,v_b,xyz_b,y_b] = get_image_data2(input_slicetime(20,:),''); % requires at least 20 volumes/files. If conversion 4D--> 3D went wrong, then error here
TR = scanTR;
nSlices = size(y_b,3);
num_slices = nSlices;

fname = '/data/home/3280411/DATA/extern/BATCHES/sub-Posval10-0060_rest.json';
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
json = jsondecode(str);
slice_times = json.SliceTiming;

[minDistance, indexOfMin] = min(abs(slice_times-(TR/2)));
reference_time = slice_times(indexOfMin);

matlabbatch{1}.spm.temporal.st.scans = {cellstr(input_slicetime)};
matlabbatch{1}.spm.temporal.st.nslices = nSlices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR - (TR/nSlices);
matlabbatch{1}.spm.temporal.st.so = [slice_times];
matlabbatch{1}.spm.temporal.st.refslice = reference_time;

% Input for realignment
input_realign_ref = [];
for ifiles = 1:size(input_slicetime,1)
    [p,f,e] = fileparts(input_slicetime(ifiles,:));
    temp1 = fullfile(p,['a',f,e]);
    input_realign_ref = strvcat(input_realign_ref,temp1);
    clear temp1 p f e
end
matlabbatch{2}.spm.spatial.realign.estwrite.data = {cellstr(input_realign_ref)};

% Input for coregistration
input_coreg_ref = [];
for ifiles = 1
    [p,f,e] = fileparts(input_realign_ref(ifiles,:));
    temp1 = fullfile(p,['mean',f,e]);
    input_coreg_ref = strvcat(input_coreg_ref,temp1);
    clear temp1 p f e
end
matlabbatch{3}.spm.spatial.coreg.estimate.ref = cellstr(input_coreg_ref);
eval(['[input_coreg_source,dirs]=spm_select(''FPList'',anatpath,''^',anatCode,'.*\.nii$'');']);
matlabbatch{3}.spm.spatial.coreg.estimate.source = cellstr(input_coreg_source);

% est & write normalise anatomy
matlabbatch{4}.spm.spatial.normalise.estwrite.subj.vol = cellstr(input_coreg_source);
matlabbatch{4}.spm.spatial.normalise.estwrite.subj.resample = cellstr(input_coreg_source);
whichSPM = which('spm');
[spmPath,~,~] = fileparts(whichSPM);
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.tpm = cellstr(fullfile(spmPath,'tpm','TPM.nii'));

% write normalise functionals using anatomy normalization
input_normfunc_def = [];
for ifiles = 1:size(input_coreg_source,1)
    [p,f,e] = fileparts(input_coreg_source(ifiles,:));
    temp1 = fullfile(p,['y_',f,e]);
    input_normfunc_def = strvcat(input_normfunc_def,temp1);
    clear temp1 p f e
end
input_normfunc_resample = [];
for ifiles = 1:size(input_realign_ref,1)
    [p,f,e] = fileparts(input_realign_ref(ifiles,:));
    temp1 = fullfile(p,['r',f,e]);
    input_normfunc_resample = strvcat(input_normfunc_resample,temp1);
    clear temp1 p f e
end

matlabbatch{5}.spm.spatial.normalise.write.subj.def = cellstr(input_normfunc_def);
matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(input_normfunc_resample);

% Input for smoothing
input_smooth = [];
for ifiles = 1:size(input_normfunc_resample,1)
    [p,f,e] = fileparts(input_normfunc_resample(ifiles,:));
    temp1 = fullfile(p,['w',f,e]);
    input_smooth = strvcat(input_smooth,temp1);
    clear temp1 p f e
end
matlabbatch{6}.spm.spatial.smooth.data = cellstr(input_smooth);

try
    spm_jobman('run_nogui',matlabbatch);
end
clear matlabbatch


%%%%% MOVEMENT ANALYSIS %%%%%
censor_scans_tmp = [];
censor_scans = [];
all_censor_scans = []; % one before the scan and two after the scan
complete_matrix = [];
matrix = [];
rpfile = dir(fullfile(funcpath,['rp_*.txt']));
rpfile = fullfile(funcpath,rpfile(1).name);

movement_matrix=load(rpfile);
movement_matrix_orig = movement_matrix; % to generate scrubbed realignment matrix for analysis with scrubbed matrix
radius = 50;
rotation_parameters =movement_matrix(:,4:6);
converted_rot_params = radius*rotation_parameters*pi/180;
movement_matrix(:,4:6) = converted_rot_params;
movement_matrix_shifted = movement_matrix(2:end,:); % beweging tussen scans
scan2scan_matrix= abs(movement_matrix(1:end-1,:)-movement_matrix_shifted); % geeft beweging
FD = sum(scan2scan_matrix,2); % sum over columns
z = ones(size(FD))'.*0.9;

% determine scans with too much movement
censor_scans_tmp = find(FD(:)>0.9);

% determine percentage of movement scans
perc_movement = (length(censor_scans_tmp)/length(movement_matrix_orig))*100;

rpmat = load(rpfile);
rpmatdiff = rpmat(2:end,:) - rpmat(1:end-1,:);
r = 65; % radius (i.e. distance between cortex and center of brain, in mm)
x = [-r 0 0 1]';
fastmotion = zeros(size(rpmatdiff,1),1);

for idyn = 1:size(rpmatdiff,1)
    q = rpmatdiff(idyn,:);
    q(4:6) = q(4:6)*pi/180;     % convert degrees into radians
    T = [1 0 0 q(1);0 1 0 q(2);0 0 1 q(3); 0 0 0 1];    % setting up the translation matrix
    % setting up the rotation matrices
    R1 = [1 0 0 0;0 cos(q(4)) sin(q(4)) 0; 0 -sin(q(4)) cos(q(4)) 0; 0 0 0 1];
    R2 = [cos(q(5)) 0 sin(q(5)) 0; 0 1 0 0; -sin(q(5)) 0 cos(q(5)) 0; 0 0 0 1];
    R3 = [cos(q(6)) sin(q(6)) 0 0;-sin(q(6)) cos(q(6)) 0 0; 0 0 1 0;0 0 0 1];
    R = R1*R2*R3;
    M = T*R;
    y = M*x;
    fastmotion(idyn) = sqrt((y(1)-x(1))^2 + (y(2)-x(2))^2 + (y(3)-x(3))^2);
    
    clear q T R* M y
end

oldnScan = length(rpmat);
%%%%% BUILDING MASKS %%%%%
%load preproc data
eval(['[ts,dirs]=spm_select(''FPList'',funcpath,''^swra',scanCode,'.*\.nii$'');']);
eval(['[t1img,dirs]=spm_select(''FPList'',anatpath,''^w',anatCode,'*'');']);
prefix = 'swra';

%The first functional...
spacedefimg=deblank(ts(1,:));

%Brain extraction
[spacedefimg]=deblank(ts(1,:));
[images,v_b,xyz_b,y_b] = get_image_data2([deblank(ts(1,:))],'message');
[p n e v] = spm_fileparts(images);
meanpath = p;
x=find(y_b(:) > (max(y_b(:))*(10/100)));

nr_voxels = size(x,1);
disp(['<mask_data> Number of voxels in mask: ',num2str(nr_voxels),'']);
mask = zeros(size(y_b));
mask(x) = ones;
v_b.pinfo(1) = 1; % to compensate for possible dcm2niix error pp
save_image_data(mask,v_b,fullfile(p,'rmask_func_file.nii'),'fu');

%Define mask variables:
maskmatrix = spm_read_vols(spm_vol(spm_select('FPList',p,['^rmask_.*.nii$'])));
mask = spm_select('FPList',p,['^rmask_.*.nii$']);

% load files:
V = spm_vol(ts);
nScan = size(V,1); %% Number of scans
xdim = V(1).dim(1); %% Dimensions ...
ydim = V(1).dim(2);
zdim = V(1).dim(3);

% Set plane coordinates
xords = (1:xdim)'*ones(1,ydim);
xords = xords(:)';
yords = ones(xdim,1)*(1:ydim);
yords = yords(:)';

%Load Mask files:
%Loads mask FILE
VM = spm_vol(mask);
% Makes a matrix from mask? again
Mask = spm_read_vols(VM);

% All the ones are the number of voxels of the mask! Why not = 1? :-)
nMVox = numel(find(Mask(:) > 0));

% Set-up high pass filter
clear HPF
clear Vhpf

% Set by me to 128...
HPF.HParam = hp_cutoff;
% Number of scans,ie episeq
HPF.row = 1:nScan;
% Time is takes for 1 scan. Total time therefore is then nScanXisi
HPF.RT = scanTR;
% Input these values to the spm_filter (High pass filter that cuts
% out all the low level noise...

HPF = spm_filter(HPF);
Vhpf(1:nScan) = deal(struct(...
    'fname',   [],...
    'dim',     V(1).dim,...
    'dt',      V(1).dt,...
    'mat',     V(1).mat,...
    'pinfo',   V(1).pinfo,...
    'descrip', ''));
for i = 1:nScan
    [p n e v] = spm_fileparts(V(i).fname);
    Vhpf(i).fname = fullfile(p,['hpf',num2str(hp_cutoff),n,e,v]);
    Vhpf(i).descrip   = [V(i).descrip,' - high-pass (',num2str(hp_cutoff),'s) filtered'];
    clear p n e v
end

%Make them into nifti's!
Vhpf = spm_create_vol(Vhpf);
[p n e v] = spm_fileparts(V(1).fname);

% check scan dimensions; important after normalization - pp
for z = 1:zdim
    % current plane-specific parameters
    zords   = z*ones(xdim*ydim,1)';
    Cm = reshape(Mask(:,:,z),[1,xdim*ydim]);
    Cm = find(Cm>0);
    if isempty(Cm) == 0
        % determine brain boundary
        if exist('lowermargin') == 0
            lowermargin = z;
        end
        uppermargin = z;
    end
end
clear Cm zords
scandim = uppermargin - lowermargin + 1;

%For timeseries graphs later on....
%The z for each scan...
nMVoxSl = zeros(1,scandim);

for i = 1:scandim
    nMVoxSl(1,i) = numel(find(Mask(:,:,i) > 0));
end

%Make Output files:
Vtmean = struct(...
    'fname',   fullfile(outputdir,['tmean_',n,e,v]),...
    'dim',     V(1).dim,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     V(1).mat,...
    'pinfo',   [1 0 0]',...
    'descrip', [V(1).descrip,' - temporal mean signal map']);
Vtmean = spm_create_vol(Vtmean);

Vtsd = struct(...
    'fname',   fullfile(outputdir,['tsd_',n,e,v]),...
    'dim',     V(1).dim,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     V(1).mat,...
    'pinfo',   [1 0 0]',...
    'descrip', [V(1).descrip,' - temporal standard deviation signal map']);
Vtsd = spm_create_vol(Vtsd);

Vtsnr = struct(...
    'fname',   fullfile(outputdir,['tsnr_',n,e,v]),...
    'dim',     V(1).dim,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     V(1).mat,...
    'pinfo',   [1 0 0]',...
    'descrip', [V(1).descrip,' - temporal signal-to-noise ratio map']);
Vtsnr = spm_create_vol(Vtsnr);

clear p n e v
clear g
clear gSF

% Grand mean scaling % pp, removed spm_global due to spikes
for i = 1:nScan
    v_b = spm_vol([V(i).fname]);
    [y_b,xyz_b] = spm_read_vols(v_b);
    y_b(y_b == 0) = NaN;
    g(i) = nanmean(y_b(:));
end
gSF = 100./g';

%Initialize ts map variables
tsrawsl = nan(nScan,scandim);
tshpfsl = nan(nScan,scandim);

clear intvals intmax integers maxintegers

%Generate Time series maps --- Signal,SD and SNR maps:
for z = 1:zdim
    % current plane-specific parameters
    zords   = z*ones(xdim*ydim,1)'; %-plane Z coordinates
    
    Cm = reshape(Mask(:,:,z),[1,xdim*ydim]);
    Cm = find(Cm>0);
    
    xyz   = [xords; yords; zords];              %-voxel coordinates
    nVox  = size(xyz,2);                        %-number of voxels in plane
    
    Y = zeros(nScan,nVox);
    for i = 1:nScan
        Y(i,Cm)  = spm_get_data(V(i).fname,xyz(:,Cm));
        if Y(i,Cm)
            sigstd(i) = std(Y(i,Cm));
            integers(i) = length(unique(Y(i,Cm)));
            maxintegers(i) = max(Y(i,Cm));
        else
            sigstd(i) = 0;
            integers(i) = 0;
            maxintegers(i) = 0;
        end
    end
    
    if exist('integers')
        intvals(z) = mean(integers);
        intmax(z) = max(maxintegers);
        intstd(z) = mean(sigstd);
    else
        intvals(z) = 0;
        intmax(z) = 0;
        intstd(z) = 0;
    end
    
    % write raw slice
    if ~isempty(Cm)
        zadj = z - lowermargin + 1;
        tsrawsl(:,zadj) = nanmean(Y(:,Cm),2);
    end
    
    for i = 1 :nScan
        jj   = NaN*ones(xdim,ydim);
        jj(:) = Y(i,:);
        spm_write_plane(Vhpf(i),jj,z);
    end
    
    % HPF
    Y = spm_filter(HPF,Y.*kron(gSF,ones(1,nVox)));
    
    % write HPF slice
    if ~isempty(Cm)
        tshpfsl(:,zadj) = nanmean(Y(:,Cm),2);
    end
    
    % Write temporal mean map
    jj   = NaN*ones(xdim,ydim);
    jj(:) = nanmean(Y,1);
    spm_write_plane(Vtmean,jj,z);
    
    % Write temporal standard deviation map
    jj   = NaN*ones(xdim,ydim);
    jj(:) = nanstd(Y,1);
    spm_write_plane(Vtsd,jj,z);
    
    % Write temporal signal-to-noise ratio map
    jj   = NaN*ones(xdim,ydim);
    jj(:) = nanmean(Y,1)./nanstd(Y,1);
    spm_write_plane(Vtsnr,jj,z);
    spm_progress_bar('Set',100*(z/zdim));
end

clear zadj Cm

% Write mean raw and hpf signal per slice per volume
tsraw = nansum(tsrawsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average raw time series across all mask voxels
tshpf = nansum(tshpfsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average hpf time series across all mask voxels

tspsc = abs(tshpf*100/nanmean(tshpf)-100);                          % Global percent signal change
tspscsl = tshpfsl.*kron(ones(nScan,1),100./nanmean(tshpfsl,1)) - 100; % Local (i.e. slice) percent signal change

[p n] = spm_fileparts(V(1).fname);
save(fullfile(outputdir,['tsrawsl_',n,'.mat']),'tsrawsl');
save(fullfile(outputdir,['tshpfsl_',n,'.mat']),'tshpfsl');
save(fullfile(outputdir,['tsraw_',n,'.mat']),'tsraw');
save(fullfile(outputdir,['tshpf_',n,'.mat']),'tshpf');
save(fullfile(outputdir,['tspsc_',n,'.mat']),'tspsc');
save(fullfile(outputdir,['tspscsl_',n,'.mat']),'tspscsl');

maptype = {'tmean' 'tsd' 'tsnr'};
for mn = 1:3
    map = dir(fullfile(outputdir,[maptype{mn},'_',n,'*.nii']));
    mapreadout = spm_vol( fullfile(outputdir,map.name) );
    mapdata = spm_read_vols( mapreadout );
    mapdata(mapdata == 0) = NaN;
    evalc([maptype{mn},'_val = nanmean(nanmean(nanmean(nanmean(mapdata))))']);
end
signal_mean = tmean_val;
signal_sd = tsd_val;
signal_snr = tsnr_val;

descriptives = sprintf(['Relative Signal Mean: ',num2str(signal_mean),'; SD: ',num2str(signal_sd),'; TSNR: ',num2str(signal_snr)]);
disp(' ')
disp(descriptives)
disp(' ')

%%%%% LOAD DATA FOR GRAPHS %%%%%
tsraw = spm_select('FPList',outputdir,['^tsraw_',n,'.*.mat$']);
load(tsraw);
tshpf = spm_select('FPList',outputdir,['^tshpf_',n,'.*.mat$']);
load(tshpf);
tspsc = spm_select('FPList',outputdir,['^tspsc_',n,'.*.mat$']);
load(tspsc);
tspscsl = spm_select('FPList',outputdir,['^tspscsl_',n,'.*.mat$']);
load(tspscsl);
cd(outputdir)

meanfile = dir(fullfile(funcpath,['mean*.nii']));
P = spm_vol(fullfile(funcpath,meanfile(1).name));
P = P(1);
fg = spm_figure('GetWin','Graphics');
fg.Position(1) = 1000;
whitebg('black')
set(gcf,'Color',[0 0 0])
set(gcf,'InvertHardcopy','off')
spm_image('Reset');
anatpos = num2cell(fg.Position);
fg.Position(4) = (cell2mat(anatpos(4))/1.25);
spm_orthviews('Image', P, [0.0 0.0 1 1]);
spm_orthviews('Xhairs','off');
figaxes = get(gcf,'children');
for y = 1:length(figaxes)
    figaxes(y).Box = 'off';
    figaxes(y).XAxis.Color = [0 0 0];
    figaxes(y).YAxis.Color = [0 0 0];
end
set(gcf,'PaperUnits','inches','Paperposition',[0 0 8 9])
i = {'10 20 30', '5 10 20', '3 5 10', '0 0 0', '-3 -5 -10', '-5 -10 -20', '-10 -20 -30'};
for sliceding = 1:length(i)
    spm_orthviews('Reposition',[str2num(i{sliceding})])
    dynjpg = strcat(['raw_dyn',num2str(sliceding),'_func.jpg']);
    print('-djpeg',dynjpg);
end
close all;

%%%%% OLDSKOOL FRONT PAGE FOR PDF WITH SIGNAL AND MOVEMENT GRAPHICS %%%%%
% truncate name if too long
nicename = func;
if length(nicename) > 30
    nicename = strcat([nicename(1:30)]);
end

fig = spm_figure('FindWin','Graphics');
if isempty(fig)
    fig = spm_figure('Create','Graphics');
end;
fig.Position(1) = 1000;
set(0,'CurrentFigure',fig);
spm_figure('Clear','Graphics');
%
%  Calculate the number of bins required
nScan = size(ts,1);
bins = 30;                                                  % Approximate maximum number of X ticks
minbin = 10;                                                % Minimum distance between X ticks
factor = ceil(nScan/(bins*minbin));                        % Factor to calculate X tick vector
% When there are many dynamics, change increments
if nScan < 50
    minbin = 5;
elseif nScan > 200
    minbin = 25;
elseif nScan > 300
    minbin = 50;
end
xtickvec = factor*minbin:factor*minbin:nScan;            % X tick vector

%Plot 1: raw and high-pass filtered time series
subplot(8,1,1)
ax = plotyy(1:nScan,tsraw,1:nScan,tshpf);

title({[datestr(now,'dd-mm-yy HH:MM'),' QC for ',nicename];},...
    'FontWeight','bold',...
    'Interpreter','none'...
    );
set(get(ax(1),'Ylabel'),'String','S (a.u.)');

set(ax(2),'YTick',[]);
set(ax(1),'XTick',xtickvec);

set(ax(1),'XLim',[0,nScan]);
set(ax(2),'XLim',[0,nScan]);

figtitle{1} = ['raw & high-pass (cutoff: ',num2str(hp_cutoff),' s) filtered signal'];
axs = get(ax);
x1 = 0.02*(axs(1).XLim(2)-axs(1).XLim(1))+axs(1).XLim(1);
y1 = 0.75*(axs(1).YLim(2)-axs(1).YLim(1))+axs(1).YLim(1);
text(x1,y1,figtitle{1},'Color','r','Interpreter','none')

%Plot 2: global percent signal change
subplot(8,1,2)
plot(tspsc)
grid on
axis tight
ylabel('dS (%)');
hgca(2) = get(gca);
figtitle{2} = ['global % change in HPF (',num2str(hp_cutoff),' s) signal'];
x1 = 0.02*(hgca(2).XLim(2)-hgca(2).XLim(1))+hgca(2).XLim(1);
y1 = 0.75*(hgca(2).YLim(2)-hgca(2).YLim(1))+hgca(2).YLim(1);
text(x1,y1,figtitle{2},'Color','r','Interpreter','none')
set(gca,'XTick',xtickvec);

% Plot 3: local (i.e. per slice) percent signal change
subplot(8,1,3:4)

clim = [0 2]; % Color bar limits

imagesc(tspscsl',clim);
colorbar('location','NorthOutside')

grid on
axis tight
ylabel('slice in dim 3');
hgca(3) = get(gca);
figtitle{3} = ['local % change in HPF (',num2str(hp_cutoff),' s) signal per slice'];
x1 = 0.02*(hgca(3).XLim(2)-hgca(3).XLim(1))+hgca(3).XLim(1);
y1 = 0.25*(hgca(3).YLim(2)-hgca(3).YLim(1))+hgca(3).YLim(1);
text(x1,y1,figtitle{3},'Color','w','Interpreter','none')
set(gca,'XTick',xtickvec);
set(gca,'XLim',[0,nScan]);

% Plot 4: realignment parameters
[rpp rpn] = spm_fileparts(deblank(ts(1,:)));

subplot(8,1,5:6)

[ax,h1,h2] = plotyy(1:oldnScan,rpmat(:,1:3),1:oldnScan,rpmat(:,4:6));
grid on

set(h1,'LineStyle','-')
set(h2,'LineStyle','--')

set(h1(1),'Color',[0 0 1]);
set(h1(2),'Color',[0 0.5 0]);
set(h1(3),'Color',[1 0 0]);

set(h2(1),'Color',[0 0 1]);
set(h2(2),'Color',[0 0.5 0]);
set(h2(3),'Color',[1 0 0]);

set(get(ax(1),'Ylabel'),'String','Translation (mm)');
set(get(ax(2),'Ylabel'),'String','Rotation (deg)');

set(ax(1),'XTick',xtickvec);
set(ax(2),'XTick',xtickvec);

set(ax(1),'XLim',[0,nScan]);
set(ax(2),'XLim',[0,nScan]);

axs = get(ax);
figtitle{4} = [rpn,' - translations & rotations'];
x1 = 0.02*(axs(1).XLim(2)-axs(1).XLim(1))+axs(1).XLim(1);
y1 = 0.75*(axs(1).YLim(2)-axs(1).YLim(1))+axs(1).YLim(1);
%text(x1,y1,figtitle{4},'Color','r','Interpreter','none');

legend('location','NorthOutside','orientation','horizontal');
[legh,objh,outh] = legend;

legh.String{1} = 'x';
legh.String{2} = 'y';
legh.String{3} = 'z';
legh.String{4} = 'pitch';
legh.String{5} = 'roll';
legh.String{6} = 'yaw';
% Plot of 'fast motion'
subplot(8,1,7:8)

plot(2:oldnScan,fastmotion);
grid on
xlabel('scan')
ylabel('motion (mm)');

hgca(5) = get(gca);
x1 = 0.02*(hgca(5).XLim(2)-hgca(5).XLim(1))+hgca(5).XLim(1);
y1 = 0.85*(hgca(5).YLim(2)-hgca(5).YLim(1))+hgca(5).YLim(1);
set(gca,'XTick',xtickvec);
set(gca,'XLim',[0,oldnScan]);

printmovement = '0';
printjunk = '0';

% Save figure
qcreport = fullfile(outputdir,['qc_report.ps']);
print('-djpeg',fullfile(outputdir,[subjectName,'_',func,'_overview.png']));
print('-append', '-dpsc2', qcreport);

clear uppermargin lowermargin

%%
%%%%% CREATE OVERLAYS FOR MEAN, SD AND SNR DATA %%%%%
% AXIAL
for imap = 1:3
    spm_figure('Clear','Graphics');
    
    SO = slover();
    
    if imap == 1
        overlay = spm_select('FPList',outputdir,['^tmean_',n,'.*\.nii$']);
    elseif imap == 2
        overlay = spm_select('FPList',outputdir,['^tsd_',n,'.*\.nii$']);
    elseif imap == 3
        overlay = spm_select('FPList',outputdir,['^tsnr_',n,'.*\.nii$']);
    end
    
    % Slice overlay settings
    SO.img(1).vol = spm_vol(t1img);                 % T1-weighted scan
    v = 2;
    SO.img(v).vol = spm_vol(overlay);               % Overlay
    
    if imap == 1
        SO.img(v).range = [0 250];                  % Color bar range of overlay
        qctype = 'mean';
    elseif imap == 2
        SO.img(v).range = [0 2];                    % Changed back to 2!
        qctype = 'sd';
    elseif imap == 3
        SO.img(v).range = [0 250];
        qctype = 'snr';
    end
    
    [cmap] = slover('getcmap', 'actc');             % Create color map
    SO.img(v).cmap = cmap;
    SO.cbar = v;                                    % Color bar represents values in map 2
    SO.transform = 'axial';                         % Slice orientation
    SO.slices = -36:12:72;                           % Coordinates of slices
    SO.printstr = 'print -noui -painters ';         % Print settings
    paint(SO);
    
    % Print JPEG for server push
    jpegname = strcat([qctype,'_ax_',nicename]);
    set(gcf,'PaperPositionMode','auto');
    print('-djpeg',jpegname);
    print('-append', '-dpsc2', qcreport);
    clear overlay SO cmap warnstr
end

%% SAGITTAL
for imap = 1:3
    spm_figure('Clear','Graphics');
    SO = slover();
    
    if imap == 1
        overlay = spm_select('FPList',outputdir,['^tmean_',n,'.*\.nii$']);
    elseif imap == 2
        overlay = spm_select('FPList',outputdir,['^tsd_',n,'.*\.nii$']);
    elseif imap == 3
        overlay = spm_select('FPList',outputdir,['^tsnr_',n,'.*\.nii$']);
    end
    
    % Slice overlay settings
    SO.img(1).vol = spm_vol(t1img);                 % T1-weighted scan
    v = 2;

    SO.img(v).vol = spm_vol(overlay);               % Overlay
    
    if imap == 1
        SO.img(v).range = [0 250];                  % Color bar range of overlay
        qctype = 'mean';
    elseif imap == 2
        SO.img(v).range = [0 2];                    % Changed back to 2!
        qctype = 'sd';
    elseif imap == 3
        SO.img(v).range = [0 250];
        qctype = 'snr';
    end
    
    [cmap] = slover('getcmap', 'actc');             % Create color map
    SO.img(v).cmap = cmap;
    SO.cbar = v;                                    % Color bar represents values in map 2
    SO.transform = 'sagittal';                      % Slice orientation
    SO.slices = -64:8:64;                           % Coordinates of slices
    SO.printstr = 'print -noui -painters ';         % Print settings
    paint(SO);
    
    jpegname = strcat([qctype,'_sag_',nicename]);
    set(gcf,'PaperPositionMode','auto');
    print('-djpeg',jpegname);
    print('-append', '-dpsc2', qcreport);
    clear overlay SO cmap warnstr hz mps;
end
evalc(['!rm ',funcpath,'/wra* ',funcpath,'/wra* ',funcpath,'/ra* ',funcpath,'/a* ',funcpath,'/hpf128*']);

close all;


