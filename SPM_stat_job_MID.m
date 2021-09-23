function SPM_stat_job_MID(analysisDir,funcDir,behavFile,TR,outputDir,sample_jobman)
warning off
%--------------------------------------------------------------------------
% Input:
%   analysisDir:    main subjectDir
%   funcDir:        directory containing preprocessed functionals
%   behavFile:      path to behavioral file
%   TR:             TR in seconds (int)
%   OutputDir:      Output for statistical results
%   sample_jobman   SPM12 template batch file
%--------------------------------------------------------------------------

% load fmri batch template
clear matlabbatch
scancodeFunc = '*'; % select all functional files in directory
load(sample_jobman);
cd(analysisDir);

% Read behav data
dat = tdfread(behavFile);
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

%% Trial: Cue 500ms | fix (6000-Cue-Star) | Star 150:450 | Resp 300 | feedback 1500

% trigger > 4 scans (80 slices, mb 2 = 40 pulses per 2.5 seconds, 160 pulses)

% create onsets per stimulus type based on target timestamp
dat.timestamp = dat.timestamp - dat.targetOnset;
dat.timestamp = dat.timestamp - dat.timestamp(expStart) + 10000; % t=0 for first trial, add 10s for dummy scans

for cue_type = 1:5 
    for row = expStart:size(dat.date,1)
        if dat.cue(row) == cue_type  & strcmp(trialcode{row},'testITI') 
            ACC{cue_type} = [ACC{cue_type} dat.ACC(row)];
            %if dat.ACC(row) == 0; continue; end; % skip incorrect
            %responses
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

cd(analysisDir);
mkdir(outputDir);
cd(outputDir);
eval('!rm SPM.mat')
analysisDir = pwd;

% define cue types and base condition labels
conType = {'anticip','outcome'};
conName = {'low_reward','high_reward','low_punishment','high_punishment','neutral'};
baseLine = {'rest','neutral'};
duration = TR; % event duration of analysis

% define functional files
eval(['[input_modelspec,dirs]=spm_select(''FPList'',funcDir,''^swra',scancodeFunc,'.*\.nii$'');']);
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(input_modelspec);
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(analysisDir);

if (length(input_modelspec)*TR) < expDuration
   warning('fMRI sequence shorter than experiment duration') % shouldn't happen of course; so maybe wrong TR entered
end

% Define base conditions
% Anticipation
offset = 2; % in s, to offset the event in relation to beginning of trial
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
iType = 2; preCond = iCond; iCond = 1;
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
i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [1 1 0 0 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward anticipation vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [1 1 0 0 -1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward anticipation vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 1 1 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment anticipation vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 1 1 -1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment anticipation vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 1 1 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward outcome vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 1 1 0 0 -1 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Reward outcome vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 0 0 1 1 0 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment outcome vs rest'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

i = i + 1;
matlabbatch{3}.spm.stats.con.consess{i}.tcon.weights = [0 0 0 0 0 0 0 1 1 -1 0];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.name = ['Punishment outcome vs neutral'];
matlabbatch{3}.spm.stats.con.consess{i}.tcon.sessrep = 'none';

% Save & run
jobfilename = ['stat_job_MID_',datestr(now,'yyyymmmdd'),'.mat'];
eval(['save ',fullfile(analysisDir,jobfilename),' matlabbatch']);
clear matlabbatch
job = fullfile(analysisDir,jobfilename);
spm_jobman('run_nogui',job);








