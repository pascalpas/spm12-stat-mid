function SPM_stat_job_MID(analysisDir,funcDir,behavDir,TR,outputDir)
warning off
%--------------------------------------------------------------------------
% Input:
%   analysisDir:    main subjectDir
%   funcDir:        directory containing preprocessed functionals
%   behavFile:      path to behavioral file
%   TR:             TR in seconds
%   OutputDir:      Output for statistical results
%--------------------------------------------------------------------------

altBehavDir = '/data/home/3280411/DATA/extern/MID_oct2022/behav' % look for files when not present in subjectdir

% load fmri batch template
clear matlabbatch
scancodeFunc = '*'; % select all functional files in directory
localPath = which('SPM_stat_job_MID.m');
[localDir,~,~] = fileparts(localPath);
sample_jobman = fullfile(localDir,'sample_jobman_stat_SPM12.mat'); 
load(sample_jobman);
cd(analysisDir);

% Read behav data
behavFile = dir(fullfile(behavDir,'*.iqdat'));
if length(behavFile) == 0
    disp('No behav found in subject directory')
    disp('>> trying to copy from other dir')
    behavFile = dir(fullfile(altBehavDir,['*raw*',analysisDir(end-6:end),'*.iqdat']));
    eval(['!cp ',fullfile(behavFile.folder,behavFile.name),' .'])
    if length(behavFile) == 0
        disp('>> No file in other dir either.. sadface')
        return;
    end
end
dateString = regexp(behavFile.name,'_202\d-\d\d-\d\d','match');
dateString = dateString{1};
dateString = dateString(2:end);
scanDate = datetime(dateString,'InputFormat','yyyy-MM-dd');
dummyDate = datetime('2021-11-25','InputFormat','yyyy-MM-dd'); % dummies changed from 4 to 8
if scanDate >= dummyDate % add offset for dummy scans
    dummyOffset = 8*TR;
else
    dummyOffset = 4*TR;
end

dat = tdfread(behavFile.name);
dat.timestamp = dat.picture0x2Etarget0x2Etimestamp; % star target; change variable name because of legibility

for cue_type = [1 2 3 4 5] % create empty variables to fill later
   onset{cue_type} = []; 
   RT{cue_type} = [];
   ACC{cue_type} = [];
end
trialcode = cellstr(dat.trialcode);
for n = 1:length(trialcode)
   if strcmp(trialcode{n},'testITI')
       expStart = n; % first trial 
       break;
   end
end
expDuration = (dat.timestamp(end) - dat.timestamp(expStart))/1000/60;

cd(analysisDir);
mkdir(outputDir);
cd(outputDir);
eval('!rm SPM.mat')
analysisDir = pwd;

% define functional files
eval(['[input_modelspec,dirs]=spm_select(''FPList'',funcDir,''^swra',scancodeFunc,'.*\.nii$'');']);
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(input_modelspec);
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(analysisDir);

% signal threshold
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8; % 0.8 default
% matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(localDir,'combined.nii'),1'}; 

if (length(input_modelspec)*TR) < (expDuration*60)
   disp('fMRI sequence shorter than experiment duration') % shouldn't happen of course; so maybe wrong TR entered
   return;
end

%% Trial: Cue 500ms | fix (6000-Cue-Star) | Star 150:450  feedback 1500
% trigger > 4 scans = 10 seconds (80 slices, mb 2 = 40 pulses per 2.5 seconds, 160 pulses)

% create onsets per stimulus type based on target timestamp
dat.timestamp = dat.timestamp - dat.targetOnset; % first timestamp is of target, so remove target onset to get actual beginning of trial in scan
dat.timestamp = dat.timestamp - dat.timestamp(expStart) + dummyOffset*1000; % t=0 for first trial, add 10s for dummy scans

for cue_type = 1:5 
    for row = expStart:size(dat.date,1)
        if dat.cue(row) == cue_type  & strcmp(trialcode{row},'testITI') 
            ACC{cue_type} = [ACC{cue_type} dat.ACC(row)];
            %if dat.ACC(row) == 0; continue; end; % uncomment to skip incorrect
            % responses
            trialOnset = dat.timestamp(row);
            try; if round(onset{cue_type}(end)) == round(trialOnset); continue; end; end % prevent duplicate entries
            onset{cue_type} = [onset{cue_type} trialOnset];   
            respRT = dat.RT(row) + dat.targetOnset(row);
            RT{cue_type} = [RT{cue_type} respRT]; % response time defined from beginning of trial
        end
    end
end
for n = 1:length(onset) 
    onsets{n} = onset{n}'./1000;  
    RT{n} = RT{n}' ./1000;
end

for cue_type = 1:5 % descriptive statistics per cue type
    behavVal(cue_type,:) = [mean(RT{cue_type})-5 mean(ACC{cue_type})];
end
%figure; bar(behavVal(:,2));

% define cue types and base condition labels
conType = {'anticip','outcome'};
conName = {'low_reward','high_reward','low_punishment','high_punishment','neutral'};
baseLine = {'rest','neutral'};

% Define base conditions
% Anticipation
offset = 0.5; % in s, to offset the event in relation to beginning of trial
duration = 5; % event duration of analysis
iType = 1; iCond = 1; 
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).name = [conType{iType} '_' conName{iCond}]; 
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).onset = onsets{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).duration = [duration];
iCond = iCond + 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).onset = onsets{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).duration = [duration];
iCond = iCond + 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).onset = onsets{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).duration = [duration];
iCond = iCond + 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).onset = onsets{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).duration = [duration];
iCond = iCond + 1; 
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).onset = onsets{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond).duration = [duration];

% Outcome
offset = 0; % in s, to offset the event in relation response
duration = TR; % event duration of analysis
iType = 2; preCond = iCond; iCond = 1;

% Remove incorrect responses from onsets
onlyCorrect = 1; % turn off when subject has very low accuracy
if onlyCorrect == 1
    for cue_type = 1:5 
        iNew = 1;
        for row = 1:length(onsets{cue_type})
            if ACC{cue_type}(row) == 1
                newOnsets{cue_type}(iNew) = onsets{cue_type}(row);
                newRT{cue_type}(iNew) = RT{cue_type}(row);
                iNew = iNew + 1;
            end
        end
    end
    onsets = newOnsets;
    RT = newRT;
end

matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).onset = onsets{iCond} + RT{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).duration = [duration];
iCond = iCond + 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).onset = onsets{iCond} + RT{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).duration = [duration];
iCond = iCond + 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).onset = onsets{iCond} + RT{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).duration = [duration];
iCond = iCond + 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).onset = onsets{iCond} + RT{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).duration = [duration];
iCond = iCond + 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).name = [conType{iType} '_' conName{iCond}];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).onset = onsets{iCond} + RT{iCond} + offset;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(iCond + preCond).duration = [duration];

% add realignment parameters to statistical analysis
realign_parameter_file =  dir(fullfile(funcDir,'rp_*.txt'));
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(fullfile(funcDir,realign_parameter_file(1).name));

% SPM.mat for model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(analysisDir,'SPM.mat'));

% SPM.mat for contrast manager
matlabbatch{3}.spm.stats.con.consess{1} = []; 
matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(analysisDir,'/SPM.mat'));

i = 1;
matlabbatch{3}.spm.stats.con.consess{i}.fcon.name = 'Effects of interest: task';
matlabbatch{3}.spm.stats.con.consess{i}.fcon.weights{1} = [ 1 0 0 0 0 0 0 0 0 0 0 0 ; ...
    0 1 0 0 0 0 0 0 0 0 0 0 ; ...
    0 0 1 0 0 0 0 0 0 0 0 0 ; ...
    0 0 0 1 0 0 0 0 0 0 0 0 ; ...
    0 0 0 0 1 0 0 0 0 0 0 0 ; ...
    0 0 0 0 0 1 0 0 0 0 0 0 ; ...
    0 0 0 0 0 0 1 0 0 0 0 0 ; ...
    0 0 0 0 0 0 0 1 0 0 0 0 ; ...
    0 0 0 0 0 0 0 0 1 0 0 0 ; ...
    0 0 0 0 0 0 0 0 0 1 0 0 ]; % matrix of all base conditions
matlabbatch{3}.spm.stats.con.consess{i}.fcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.fcon.name = 'Effect of motion'; % 10 base conditions plus 6 realignment parameters
matlabbatch{3}.spm.stats.con.consess{i}.fcon.weights{1} = [0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0];
matlabbatch{3}.spm.stats.con.consess{i}.fcon.sessrep = 'none';

% tcons | conditions: low_reward, high_reward, low_punishment, high_punishment, neutral
i = i + 1; % 0003
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0.5 0.5 0 0 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward anticipation vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0004
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0.5 0.5 0 0 -1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward anticipation vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0005
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0.5 0.5 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment anticipation vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0006
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0.5 0.5 -1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment anticipation vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0007
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 0.5 0.5 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward outcome vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0008
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 0.5 0.5 0 0 -1 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward outcome vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0009
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 0 0 0.5 0.5 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment outcome vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0010
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 0 0 0.5 0.5 -1 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment outcome vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0011
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [1 1 1 1 1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Anticipation all vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1; % 0012
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 1 1 1 1 1 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Outcome all vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

% Save & run
jobfilename = ['stat_job_MID_',datestr(now,'yyyymmmdd'),'.mat'];
eval(['save ',fullfile(analysisDir,jobfilename),' matlabbatch']);
clear matlabbatch
job = fullfile(analysisDir,jobfilename);
spm_jobman('run_nogui',job);








