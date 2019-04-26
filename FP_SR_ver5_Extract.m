% Trim function
% plot - subplot 20181221
% Manual Trial inserted

%% Extract photometry data from TDT Data Tank

%path_to_dropbox = '/Users/marcfuccillo/Desktop/Opeyemi Photometry Analysis'; % point to folder where data is
%tankdir = path_to_dropbox;
%tankname = 'SSC312-171018-134706';
%blockname = 'SSC312_1024';
%new_direc = strcat(tankdir, '\',tankname, '\', blockname);
%cd (new_direc);
%storenames = {'470B' '405B' 'Ms1/' 'Ms1\' 'Ms2/' 'Ms2\' };



close all;
clear all;

Trim_Start = 0; % Trim window at the first part 
Trim_End = 30000; % Trim window at the last part


filepath=uigetdir;
%exptname=base; filename=dir+exptname; savename=filename+100Hz
filepath_sep = regexp(filepath, filesep, 'split');
expnameDef = filepath_sep(end);
exptname = filepath_sep(end);
exptname = exptname{1};
%exptname = inputdlg('experiment name?','Input',1,filepath_sep(end),'on'); exptname = exptname{1};
mkdir(exptname); addpath(exptname);% make a folder name for the m file output
%     exptnameDef=strsplit(filepath,'/');
%     exptnameDef=char(exptnameDef(length(exptnameDef)));  %the savename becomes the folder name
%     exptnameDef = {exptnameDef}; %default expt name saved as a cell array of strings
%     filename=strcat(exptname,'/',exptname); %savename dir becomes destination and filename


filename = strcat(exptname,'.mat'); %filename where entire workspace can be saved
figname = strcat(exptname,'.fig');  %figname for three signal graphs
figname1 = strcat(exptname,'_final','.fig');  %figname for final data graph
storenames = {'470P' '405P' 'Tal/' 'Rrd/' '41P2' '46P2'}; % '41P2' '46P2'
%
% Channel_List = {'Channel_1','Channel_2'};
% [ChannelLookup,tf] = listdlg('ListString',Channel_List);

% if ChannelLookup == 1
%     storenames = {'470P' '405P' 'Tal/' 'Rrd/' }; % '41P2' '46P2'
% elseif ChannelLookup == 2
%     storenames = {'41P2' '46P2' 'Tal/' 'Rrd/' };
% end

uiopen

%% extract

for k = 1:numel(storenames)
    storename = storenames{k};
    S{k} = TDT_MATconvertBDH(filepath, storename); %BHedit
end

%% Massage data and get time stamps
% Get LMag data as a vector (repeat for each channel)
dat1_ch1 = S{1}; %470nm data
dat1_ch1.data = reshape(dat1_ch1.data', [],1); % unwrap data from m x 256 array
dat1_ch1.data = dat1_ch1.data(1+Trim_Start:end-Trim_End); %trim off end, why?
dat2_ch1 = S{2}; %405nm data
dat2_ch1.data = reshape(dat2_ch1.data', [],1); % unwrap data from m x 256 array
dat2_ch1.data = dat2_ch1.data(1+Trim_Start:end-Trim_End);

dat1_ch2 = S{5}; %470nm data
dat1_ch2.data = reshape(dat1_ch2.data', [],1); % unwrap data from m x 256 array
dat1_ch2.data = dat1_ch2.data(1+Trim_Start:end-Trim_End); %trim off end, why?
dat2_ch2 = S{6}; %405nm data
dat2_ch2.data = reshape(dat2_ch2.data', [],1); % unwrap data from m x 256 array
dat2_ch2.data = dat2_ch2.data(1+Trim_Start:end-Trim_End);



% Get LMag timestamps (use chani1 - timestamps should be the same for all LMag channels
ts = dat1_ch1.timestamps;
t_rec_start = ts(1);
samplingRate = dat1_ch1.sampling_rate;
samplingRate_ch1 = dat1_ch1.sampling_rate; % TDT returns 'true' sample rate.

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:dat1_ch1.npoints-1)*(1./dat1_ch1.sampling_rate));
ts = reshape(ts',[],1);
ts = ts(1+Trim_Start:end-Trim_End);

% % Get TTL input timestamps
CorrectPress = S{4}.timestamps - t_rec_start;
% Initiation = S{3}.timestamps - t_rec_start;
Trial_LitOn = S{3}.timestamps - t_rec_start;
% RewardRetrieval = S{5}.timestamps - t_rec_start;

%% Single fit - if both channels change consistently for the whole session

% Smooth dat1 and dat2, fit dat2 (control) to dat1 (signal)

dat1_ch1 = filtfilt(ones(1,100)/100,1, dat1_ch1.data); %MF:digital filter
dat2_ch1 = filtfilt(ones(1,100)/100,1, dat2_ch1.data);

tsOrig=ts; %MF:original tstamps
[dat1_ch1,ts]=resample(dat1_ch1,tsOrig,40,10,10); %MF:downsampling from ~1000Hz to 40Hz, ts=new SR
dat2_ch1=resample(dat2_ch1,tsOrig,40,10,10);
samplingRate = 40;

debleach_flag=1; %MF:decays to 0
[dfof, mod_dat1_ch1, fit_dat1_ch1, offset_dat1_ch1]=debleachBDH(ts',dat1_ch1,debleach_flag);
[dfof_control, mod_dat2_ch1, fit_dat2_ch1, offset_dat2_ch1] =debleachBDH(ts',dat2_ch1,debleach_flag);
%MF: this dfof_control is a problem


%% Channel1
dfofCorr=subtract_refBDH(ts',dfof,dfof_control,'None'); %MF:
dfofCorr_sub=subtract_refBDH(ts',dfof,dfof_control,'Subtract');

[controlFit_ch1] = controlFit_ch1 (dat1_ch1, dat2_ch1);
% from here on, what ope did was take three different ways to calculate dfof,
% there was the one from ES originally, that is the normdat from her code.
% Anything with just normdat with no db is that. anything with db is debleached.
% z means it got z scored, with no moving window. modz means moving window.
% Anythign with diff is a differenced time series.
% Get delta F/F using controlFit
[normDat_ch1] = deltaFF (dat1_ch1, controlFit_ch1); %MF:controlFit is weird, making [normDat] weird

%global z-score with 470ch only
normDat_DB_ch1 = 100*dfofCorr'; %MF:looks OK, doesn't use 405 for subtraction
normDat_DBZ_ch1 = (normDat_DB_ch1-mean(normDat_DB_ch1))/std(normDat_DB_ch1); %MF:OK, potential option.
normDat_DBZ_diff_ch1 = diff(normDat_DBZ_ch1); %MF: forget all difference measures


%global z-score with 470ch-405ch, problem is 405ch is much larger dfof;
normDat_DB_sub_ch1 = 100*dfofCorr_sub'; %MF: can't use this, 405 overwhelms
normDat_DBZ_sub_ch1 = (normDat_DB_sub_ch1-mean(normDat_DB_sub_ch1))/std(normDat_DB_sub_ch1); %MF: nope
normDat_DBZ_sub_diff_ch1 = diff(normDat_DBZ_sub_ch1); %MF:nope

%calculate local means and std in 2 minute windows, MF: all done on
%[normDat], which is not ideal
moving_mean_ch1 = zeros(1,length(normDat_ch1));
moving_sd_ch1 = zeros(1,length(normDat_ch1));
moving_mean_ch1(1:120*samplingRate) = mean(normDat_ch1(1:120*samplingRate));
moving_sd_ch1(1:120*samplingRate) = std(normDat_ch1(1:120*samplingRate));
moving_mean_ch1(end-(120*samplingRate):end) = mean(normDat_ch1(end-(120*samplingRate):end));
moving_sd_ch1(end-(120*samplingRate):end) = std(normDat_ch1(end-(120*samplingRate):end));
for i = (120*samplingRate)+1:(length(normDat_ch1)-(120*samplingRate))
    moving_mean_ch1(i) = mean(normDat_ch1((i-120*samplingRate):(i+120*samplingRate)));
    moving_sd_ch1(i) = std(normDat_ch1((i-120*samplingRate):(i+120*samplingRate)));
end

%z-z-score data using local means and standard deviations
% MF: mix of normDat and normDat_DB
normDat_modZ_ch1 = zeros(1,length(normDat_DB_ch1));
for i = 1:length(normDat_ch1)
    normDat_modZ_ch1(i) = (normDat_ch1(i)-moving_mean_ch1(i))/moving_sd_ch1(i);
end

normDat_modZ_ch1 = normDat_modZ_ch1';

%calculate local means and std in 2 minute windows for debleached data
%with no subtraction
moving_mean_ch1 = zeros(1,length(normDat_DB_ch1));
moving_sd_ch1 = zeros(1,length(normDat_DB_ch1));
moving_mean_ch1(1:120*samplingRate) = mean(normDat_DB_ch1(1:120*samplingRate));
moving_sd_ch1(1:120*samplingRate) = std(normDat_DB_ch1(1:120*samplingRate));
moving_mean_ch1(end-(120*samplingRate):end) = mean(normDat_DB_ch1(end-(120*samplingRate):end));
moving_sd_ch1(end-(120*samplingRate):end) = std(normDat_DB_ch1(end-(120*samplingRate):end));
for i = (120*samplingRate)+1:(length(normDat_DB_ch1)-(120*samplingRate))
    moving_mean_ch1(i) = mean(normDat_DB_ch1((i-120*samplingRate):(i+120*samplingRate)));
    moving_sd_ch1(i) = std(normDat_DB_ch1((i-120*samplingRate):(i+120*samplingRate)));
end

%z-score data using local means and standard deviations, just 470ch
normDat_DB_modZ_ch1 = zeros(1,length(normDat_DB_ch1));
for i = 1:length(normDat_DB_ch1)
    normDat_DB_modZ_ch1(i) = (normDat_DB_ch1(i)-moving_mean_ch1(i))/moving_sd_ch1(i);
end

normDat_DB_modZ_ch1 = normDat_DB_modZ_ch1';

%calculate local means and std in 2 minute windows for debleached data with
%subtraction of 405
moving_mean_ch1 = zeros(1,length(normDat_DB_sub_ch1));
moving_sd_ch1 = zeros(1,length(normDat_DB_sub_ch1));
moving_mean_ch1(1:120*samplingRate) = mean(normDat_DB_sub_ch1(1:120*samplingRate));
moving_sd_ch1(1:120*samplingRate) = std(normDat_DB_sub_ch1(1:120*samplingRate));
moving_mean_ch1(end-(120*samplingRate):end) = mean(normDat_DB_sub_ch1(end-(120*samplingRate):end));
moving_sd_ch1(end-(120*samplingRate):end) = std(normDat_DB_sub_ch1(end-(120*samplingRate):end));
for i = (120*samplingRate)+1:(length(normDat_DB_sub_ch1)-(120*samplingRate))
    moving_mean_ch1(i) = mean(normDat_DB_sub_ch1((i-120*samplingRate):(i+120*samplingRate)));
    moving_sd_ch1(i) = std(normDat_DB_sub_ch1((i-120*samplingRate):(i+120*samplingRate)));
end

%z-score data using local means and standard deviations
normDat_DB_sub_modZ_ch1 = zeros(1,length(normDat_DB_sub_ch1));
for i = 1:length(normDat_DB_sub_ch1)
    normDat_DB_sub_modZ_ch1(i) = (normDat_DB_sub_ch1(i)-moving_mean_ch1(i))/moving_sd_ch1(i);
end

normDat_DB_sub_modZ_ch1 = normDat_DB_sub_modZ_ch1';



%% Massage data and get time stamps
% Get LMag data as a vector (repeat for each channel)
dat1_ch2 = S{5}; %470nm data
dat1_ch2.data = reshape(dat1_ch2.data', [],1); % unwrap data from m x 256 array
dat1_ch2.data = dat1_ch2.data(1+Trim_Start:end-Trim_End); %trim off end, why?
dat2_ch2 = S{6}; %405nm data
dat2_ch2.data = reshape(dat2_ch2.data', [],1); % unwrap data from m x 256 array
dat2_ch2.data = dat2_ch2.data(1+Trim_Start:end-Trim_End);

% Get LMag timestamps (use chani1 - timestamps should be the same for all LMag channels
ts = dat1_ch2.timestamps;
t_rec_start = ts(1);
samplingRate = dat1_ch2.sampling_rate;
samplingRate_ch2 = dat1_ch2.sampling_rate; % TDT returns 'true' sample rate.

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:dat1_ch2.npoints-1)*(1./dat1_ch2.sampling_rate));
ts = reshape(ts',[],1);
ts = ts(1+Trim_Start:end-Trim_End);

% % Get TTL input timestamps
CorrectPress = S{4}.timestamps - t_rec_start;
% Initiation = S{3}.timestamps - t_rec_start;
Trial_LitOn = S{3}.timestamps - t_rec_start;
% RewardRetrieval = S{5}.timestamps - t_rec_start;

%% Single fit - if both channels change consistently for the whole session

% Smooth dat1 and dat2, fit dat2 (control) to dat1 (signal)

dat1_ch2 = filtfilt(ones(1,100)/100,1, dat1_ch2.data); %MF:digital filter
dat2_ch2 = filtfilt(ones(1,100)/100,1, dat2_ch2.data);

tsOrig=ts; %MF:original tstamps
[dat1_ch2,ts]=resample(dat1_ch2,tsOrig,40,10,10); %MF:downsampling from ~1000Hz to 40Hz, ts=new SR
dat2_ch2=resample(dat2_ch2,tsOrig,40,10,10);
samplingRate = 40;

debleach_flag=1; %MF:decays to 0
[dfof, mod_dat1_ch2, fit_dat1_ch2, offset_dat1_ch2]=debleachBDH(ts',dat1_ch2,debleach_flag);
[dfof_control, mod_dat2_ch2, fit_dat2_ch2, offset_dat2_ch2] =debleachBDH(ts',dat2_ch2,debleach_flag);
%MF: this dfof_control is a problem




%%
dfofCorr=subtract_refBDH(ts',dfof,dfof_control,'None'); %MF:
dfofCorr_sub=subtract_refBDH(ts',dfof,dfof_control,'Subtract');

[controlFit_ch2] = controlFit_ch2 (dat1_ch2, dat2_ch2);
% from here on, what ope did was take three different ways to calculate dfof,
% there was the one from ES originally, that is the normdat from her code.
% Anything with just normdat with no db is that. anything with db is debleached.
% z means it got z scored, with no moving window. modz means moving window.
% Anythign with diff is a differenced time series.
% Get delta F/F using controlFit
[normDat_ch2] = deltaFF (dat1_ch2, controlFit_ch2); %MF:controlFit is weird, making [normDat] weird

%global z-score with 470ch only
normDat_DB_ch2 = 100*dfofCorr'; %MF:looks OK, doesn't use 405 for subtraction
normDat_DBZ_ch2 = (normDat_DB_ch2-mean(normDat_DB_ch2))/std(normDat_DB_ch2); %MF:OK, potential option.
normDat_DBZ_diff_ch2 = diff(normDat_DBZ_ch2); %MF: forget all difference measures

%global z-score with 470ch-405ch, problem is 405ch is much larger dfof;
normDat_DB_sub_ch2 = 100*dfofCorr_sub'; %MF: can't use this, 405 overwhelms
normDat_DBZ_sub_ch2 = (normDat_DB_sub_ch2-mean(normDat_DB_sub_ch2))/std(normDat_DB_sub_ch2); %MF: nope
normDat_DBZ_sub_diff_ch2 = diff(normDat_DBZ_sub_ch2); %MF:nope

%calculate local means and std in 2 minute windows, MF: all done on
%[normDat], which is not ideal
moving_mean_ch2 = zeros(1,length(normDat_ch2));
moving_sd_ch2 = zeros(1,length(normDat_ch2));
moving_mean_ch2(1:120*samplingRate) = mean(normDat_ch2(1:120*samplingRate));
moving_sd_ch2(1:120*samplingRate) = std(normDat_ch2(1:120*samplingRate));
moving_mean_ch2(end-(120*samplingRate):end) = mean(normDat_ch2(end-(120*samplingRate):end));
moving_sd_ch2(end-(120*samplingRate):end) = std(normDat_ch2(end-(120*samplingRate):end));
for i = (120*samplingRate)+1:(length(normDat_ch2)-(120*samplingRate))
    moving_mean_ch2(i) = mean(normDat_ch2((i-120*samplingRate):(i+120*samplingRate)));
    moving_sd_ch2(i) = std(normDat_ch2((i-120*samplingRate):(i+120*samplingRate)));
end

%z-z-score data using local means and standard deviations
% MF: mix of normDat and normDat_DB
normDat_modZ_ch2 = zeros(1,length(normDat_DB_ch2));
for i = 1:length(normDat_ch2)
    normDat_modZ_ch2(i) = (normDat_ch2(i)-moving_mean_ch2(i))/moving_sd_ch2(i);
end

normDat_modZ_ch2 = normDat_modZ_ch2';

%calculate local means and std in 2 minute windows for debleached data
%with no subtraction
moving_mean_ch2 = zeros(1,length(normDat_DB_ch2));
moving_sd_ch2 = zeros(1,length(normDat_DB_ch2));
moving_mean_ch2(1:120*samplingRate) = mean(normDat_DB_ch2(1:120*samplingRate));
moving_sd_ch2(1:120*samplingRate) = std(normDat_DB_ch2(1:120*samplingRate));
moving_mean_ch2(end-(120*samplingRate):end) = mean(normDat_DB_ch2(end-(120*samplingRate):end));
moving_sd_ch2(end-(120*samplingRate):end) = std(normDat_DB_ch2(end-(120*samplingRate):end));
for i = (120*samplingRate)+1:(length(normDat_DB_ch2)-(120*samplingRate))
    moving_mean_ch2(i) = mean(normDat_DB_ch2((i-120*samplingRate):(i+120*samplingRate)));
    moving_sd_ch2(i) = std(normDat_DB_ch2((i-120*samplingRate):(i+120*samplingRate)));
end

%z-score data using local means and standard deviations, just 470ch
normDat_DB_modZ_ch2 = zeros(1,length(normDat_DB_ch2));
for i = 1:length(normDat_DB_ch2)
    normDat_DB_modZ_ch2(i) = (normDat_DB_ch2(i)-moving_mean_ch2(i))/moving_sd_ch2(i);
end

normDat_DB_modZ_ch2 = normDat_DB_modZ_ch2';

%calculate local means and std in 2 minute windows for debleached data with
%subtraction of 405
moving_mean_ch2 = zeros(1,length(normDat_DB_sub_ch2));
moving_sd_ch2 = zeros(1,length(normDat_DB_sub_ch2));
moving_mean_ch2(1:120*samplingRate) = mean(normDat_DB_sub_ch2(1:120*samplingRate));
moving_sd_ch2(1:120*samplingRate) = std(normDat_DB_sub_ch2(1:120*samplingRate));
moving_mean_ch2(end-(120*samplingRate):end) = mean(normDat_DB_sub_ch2(end-(120*samplingRate):end));
moving_sd_ch2(end-(120*samplingRate):end) = std(normDat_DB_sub_ch2(end-(120*samplingRate):end));
for i = (120*samplingRate)+1:(length(normDat_DB_sub_ch2)-(120*samplingRate))
    moving_mean_ch2(i) = mean(normDat_DB_sub_ch2((i-120*samplingRate):(i+120*samplingRate)));
    moving_sd_ch2(i) = std(normDat_DB_sub_ch2((i-120*samplingRate):(i+120*samplingRate)));
end

%z-score data using local means and standard deviations
normDat_DB_sub_modZ_ch2 = zeros(1,length(normDat_DB_sub_ch2));
for i = 1:length(normDat_DB_sub_ch2)
    normDat_DB_sub_modZ_ch2(i) = (normDat_DB_sub_ch2(i)-moving_mean_ch2(i))/moving_sd_ch2(i);
end

normDat_DB_sub_modZ_ch2 = normDat_DB_sub_modZ_ch2';



%%
f(1) = figure('Color',[1 1 1]);
% zoom on;
subplot (2, 1, 1);
plot(ts/60,dat1_ch1);
hold on
plot(ts/60,fit_dat1_ch1+offset_dat1_ch1,'r');
plot(ts/60,mod_dat1_ch1,'g');
legend('Signal','Exponential Fit','Debleached');
legend BOXOFF;
xlabel('Time(m)');
ylabel('Fluorescence');
title('Ch1 Debleaching Results for 470');
% f(2) = figure('Color',[1 1 1]);
% zoom on;
subplot (2, 1, 2);
plot(ts/60,dat2_ch1);
hold on
plot(ts/60,fit_dat2_ch1+offset_dat2_ch1,'r');
plot(ts/60,mod_dat2_ch1,'g');
legend('Signal','Exponential Fit','Debleached');
legend BOXOFF;
xlabel('Time(m)');
ylabel('Fluorescence');
title('Ch1 Debleaching Results for 405');

%
f(2) = figure('Color',[1 1 1]);
% zoom on;
subplot (2, 1, 1);
plot(ts/60,dat1_ch2);
hold on
plot(ts/60,fit_dat1_ch2+offset_dat1_ch2,'r');
plot(ts/60,mod_dat1_ch2,'g');
legend('Signal','Exponential Fit','Debleached');
legend BOXOFF;
xlabel('Time(m)');
ylabel('Fluorescence');
title('Ch2 Debleaching Results for 470');

% f(2) = figure('Color',[1 1 1]);
% zoom on;
subplot (2, 1, 2);
plot(ts/60,dat2_ch2);
hold on
plot(ts/60,fit_dat2_ch2+offset_dat2_ch2,'r');
plot(ts/60,mod_dat2_ch2,'g');
legend('Signal','Exponential Fit','Debleached');
legend BOXOFF;
xlabel('Time(s)');
ylabel('Fluorescence');
title('Ch2 Debleaching Results for 405');

%% Channel 1 review Plot
%%%%%plotting of all versions happens below
%this is old. it will plot everything, thats just so that you can inspect
%it to see what looks good. what you can tell a lot is whether the
%subtracting 405 adds noise, whether you want to, etc.
f(3) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 1);
plot(ts/60, normDat_ch1); %normDat data - not good
xlabel('Time(m)');
ylabel('\DeltaF/F (%)');
title('Original'); %No Transformations
%
% f(4) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 2);
plot(ts/60, normDat_modZ_ch1); %uses normDat data - not good
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Local Mean/SD); Original)');

% f(5) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 3);
plot(ts/60, normDat_DB_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (%)');
title('normDat_DB');
%
% f(6) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 4);
plot(ts/60, normDat_DBZ_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Global Mean/SD; Debleached)');

% f(7) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 5);
plot(ts/60, normDat_DB_modZ_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Local Mean/SD); Debleached)');

% f(8) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 6);
plot(ts(2:end)/60, normDat_DBZ_diff_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Diff Debleached Z-Score)');

% f(9) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 7);
plot(ts/60, normDat_DB_sub_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (%)');
title('Debleached Subt 405');

% f(10) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 8);
plot(ts/60, normDat_DBZ_sub_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Z-S (Global SD); Debleached w/ 405 Subt');

% f(11) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 9);
plot(ts/60, normDat_DB_sub_modZ_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Z-Score (Local SD); Debleached w/ 405 Subt');
%
% f(12) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 10);
plot(ts(2:end)/60, normDat_DBZ_sub_diff_ch1);
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Diff Debleached Z-Score w/ 405 Subt');

%% Channel 2 review plot
%%%%%plotting of all versions happens below
%this is old. it will plot everything, thats just so that you can inspect
%it to see what looks good. what you can tell a lot is whether the
%subtracting 405 adds noise, whether you want to, etc.
f(4) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 1);
plot(ts/60, normDat_ch2); %normDat data - not good
xlabel('Time(m)');
ylabel('\DeltaF/F (%)');
title('Original'); %No Transformations
%
% f(4) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 2);
plot(ts/60, normDat_modZ_ch2); %uses normDat data - not good
xlabel('Time(m)');
ylabel('\DeltaF/F (z-score)');
title('Local Mean/SD); Original)');

% f(5) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 3);
plot(ts/60, normDat_DB_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (%)');
title('normDat_DB');
%
% f(6) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 4);
plot(ts/60, normDat_DBZ_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (z-score)');
title('Global Mean/SD; Debleached)');

% f(7) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 5);
plot(ts/60, normDat_DB_modZ_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (z-score)');
title('Local Mean/SD); Debleached)');

% f(8) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 6);
plot(ts(2:end)/60, normDat_DBZ_diff_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (z-score)');
title('Diff Debleached Z-Score)');

% f(9) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 7);
plot(ts/60, normDat_DB_sub_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (%)');
title('Debleached Subt 405');

% f(10) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 8);
plot(ts/60, normDat_DBZ_sub_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (z-score)');
title('Z-S (Global SD); Debleached w/ 405 Subt');

% f(11) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 9);
plot(ts/60, normDat_DB_sub_modZ_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (z-score)');
title('Z-Score (Local SD); Debleached w/ 405 Subt');
%
% f(12) = figure('Color',[1 1 1]);
zoom on;
subplot (5, 2, 10);
plot(ts(2:end)/60, normDat_DBZ_sub_diff_ch2);
xlabel('Time(s)');
ylabel('\DeltaF/F (z-score)');
title('Diff Debleached Z-Score w/ 405 Subt');
%% Save raw, fit and dFF figs to single file

% savefig (gcf, figname);

%% Bpod loading & tagging

SessionData.Omission = cell2mat(SessionData.Omission);


%% Reward_retrieval including unrewarded trial (nan in case there is no reward retrieval)
for currentTrial = 1:SessionData.nTrials
%     [m,n] = size(SessionData.RawEvents.Trial{currentTrial}.States.Initiation_Sustain);
%     raw_scan{currentTrial} = m * n;
%     %raw_scan{currentTrial} = length(SessionData.RawEvents.Trial{currentTrial}.States.Initiation_Sustain)

        for nEvents = 1:length(SessionData.RawEvents.Trial{currentTrial}.Events.Port2In)
%             extract_scan{nEvents} = NaN(length(SessionData.RawEvents.Trial{currentTrial}.Events.Port2In));
            extract_scan = SessionData.RawEvents.Trial{currentTrial}.Events.Port2In(nEvents)-SessionData.RawEvents.Trial{currentTrial}.States.Choice(end,end);
            if extract_scan > 0
                SessionData.RawEvents.Trial{currentTrial}.States.Reward_Retrieval = SessionData.RawEvents.Trial{currentTrial}.Events.Port2In(nEvents);
                break
            else
                SessionData.RawEvents.Trial{currentTrial}.States.Reward_Retrieval = nan;
            end
            
        end
       
        Reward_Retrieval(currentTrial) = SessionData.RawEvents.Trial{currentTrial}.States.Reward_Retrieval;
end


% Block End/Start tag
Trial_BlockEnd = cumsum(SessionData.BlockLength);
Trial_BlockStart_1 = Trial_BlockEnd + 1;
Trial_BlockStart_2 = Trial_BlockEnd + 2;


SessionData.Tag.Trial_BlockEnd = zeros(1, SessionData.nTrials);
SessionData.Tag.Trial_BlockStart_1 = zeros(1, SessionData.nTrials);
SessionData.Tag.Trial_BlockStart_2 = zeros(1, SessionData.nTrials);
SessionData.Tag.Trial_BlockInflex = zeros(1, SessionData.nTrials);

for currentBlock = 1:length(Trial_BlockEnd)
    SessionData.Tag.Trial_BlockEnd(1, Trial_BlockEnd(currentBlock)) = 1;
    SessionData.Tag.Trial_BlockStart_1(1, Trial_BlockStart_1(currentBlock)) = 1;
    SessionData.Tag.Trial_BlockStart_2(1, Trial_BlockStart_2(currentBlock)) = 1;
end

Manual_Trial = []; %pre-assign


%% Extract Latency data from Bpod State
% Trial sorting by results

%TimeStamp from Bpod (20181128)

for currentTrial= 1:SessionData.nTrials
    SessionData.TimeStamp.CenterLitOn(currentTrial) = SessionData.RawEvents.Trial{currentTrial}.States.Initiation(end, 1);
    SessionData.TimeStamp.Initiation(currentTrial) = SessionData.RawEvents.Trial{currentTrial}.States.Initiation(end, 2);
    SessionData.TimeStamp.Initiation_Duration(currentTrial) = SessionData.RawEvents.Trial{currentTrial}.States.Initiation_Sustain(end,2)- SessionData.RawEvents.Trial{currentTrial}.States.Initiation_Sustain(1,1);
    SessionData.TimeStamp.Choice(currentTrial) = SessionData.RawEvents.Trial{currentTrial}.States.Choice(end, 2);
    SessionData.TimeStamp.RewardRetrieval(currentTrial) = SessionData.RawEvents.Trial{currentTrial}.States.Reward_Retrieval(1);
    SessionData.TimeStamp.ISI(currentTrial) = 0;
end

% FiberPhotometry TimeStamp Calculation from Bpod : Confirm
% TimeStamp_FP.CorrectPress_raw = CorrectPress' - (Trim_Start/samplingRate_ch1);
% for currentTrial= 1:SessionData.nTrials 
%     TimeStamp_FP.CenterLitOn(currentTrial) = TimeStamp_FP.CorrectPress_raw(currentTrial) - SessionData.TimeStamp.Choice(currentTrial) + SessionData.TimeStamp.CenterLitOn(currentTrial);
%     TimeStamp_FP.Initiation(currentTrial) = TimeStamp_FP.CorrectPress_raw(currentTrial) - SessionData.TimeStamp.Choice(currentTrial) + SessionData.TimeStamp.Initiation(currentTrial);
%     TimeStamp_FP.Choice(currentTrial) = TimeStamp_FP.CorrectPress_raw(currentTrial);
%     TimeStamp_FP.RewardRetrieval(currentTrial) = TimeStamp_FP.CorrectPress_raw(currentTrial) - SessionData.TimeStamp.Choice(currentTrial) + SessionData.TimeStamp.RewardRetrieval(currentTrial);
%     TimeStamp_FP.ISI(currentTrial) = TimeStamp_FP.CorrectPress_raw(currentTrial) - SessionData.TimeStamp.Choice(currentTrial) + SessionData.TimeStamp.ISI(currentTrial);
% end

% FiberPhotometry TimeStamp Calculation from Bpod : from Center LitOn
% TimeStamp_FP.Trial_LitOn_raw = Trial_LitOn' - (Trim_Start/samplingRate_ch1);
% for currentTrial= 1:SessionData.nTrials 
%     TimeStamp_FP.CenterLitOn(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial) ; 
%     TimeStamp_FP.Initiation(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial) + SessionData.TimeStamp.Initiation(currentTrial) - SessionData.TimeStamp.CenterLitOn(currentTrial);
%     TimeStamp_FP.Choice(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial) + SessionData.TimeStamp.Choice(currentTrial)- SessionData.TimeStamp.CenterLitOn(currentTrial);
%     TimeStamp_FP.RewardRetrieval(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial) + SessionData.TimeStamp.RewardRetrieval(currentTrial)- SessionData.TimeStamp.CenterLitOn(currentTrial);
%     TimeStamp_FP.ISI(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial) + SessionData.TimeStamp.RewardRetrieval(currentTrial) + 5- SessionData.TimeStamp.CenterLitOn(currentTrial);
% end

TimeStamp_FP.Trial_LitOn_raw = Trial_LitOn' - (Trim_Start/samplingRate_ch1);
for currentTrial= 1:SessionData.nTrials 
    Omission_Counter = sum(SessionData.Omission(1:currentTrial));
    TimeStamp_FP.CenterLitOn(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial + Omission_Counter) ; 
    TimeStamp_FP.Initiation(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial + Omission_Counter) + SessionData.TimeStamp.Initiation(currentTrial) - SessionData.TimeStamp.CenterLitOn(currentTrial);
    TimeStamp_FP.Choice(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial + Omission_Counter) + SessionData.TimeStamp.Choice(currentTrial)- SessionData.TimeStamp.CenterLitOn(currentTrial);
    TimeStamp_FP.RewardRetrieval(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial + Omission_Counter) + SessionData.TimeStamp.RewardRetrieval(currentTrial)- SessionData.TimeStamp.CenterLitOn(currentTrial);
    TimeStamp_FP.ISI(currentTrial) = TimeStamp_FP.Trial_LitOn_raw(currentTrial + Omission_Counter) + SessionData.TimeStamp.RewardRetrieval(currentTrial) + 5- SessionData.TimeStamp.CenterLitOn(currentTrial);
end


%% All Trial

for currentTrial= 1:SessionData.nTrials
    TimeStamp_AllTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
    TimeStamp_AllTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
    TimeStamp_AllTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
    TimeStamp_AllTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    
end

TimeStamp_AllTrial.CenterLitOn(isnan(TimeStamp_AllTrial.CenterLitOn))=[];
TimeStamp_AllTrial.Initiation(isnan(TimeStamp_AllTrial.Initiation))=[];
TimeStamp_AllTrial.Choice(isnan(TimeStamp_AllTrial.Choice))=[];
TimeStamp_AllTrial.RewardRetrieval(isnan(TimeStamp_AllTrial.RewardRetrieval))=[];


% Rewarded Trial

for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == 1
        Reward_Trial(currentTrial) = currentTrial;
    else
        Reward_Trial(currentTrial) = nan;
    end
    
end
Reward_Trial(isnan(Reward_Trial))=[];




for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == 1
        TimeStamp_RewardedTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_RewardedTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_RewardedTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_RewardedTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_RewardedTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardedTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardedTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardedTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_RewardedTrial.CenterLitOn(isnan(TimeStamp_RewardedTrial.CenterLitOn))=[];
TimeStamp_RewardedTrial.Initiation(isnan(TimeStamp_RewardedTrial.Initiation))=[];
TimeStamp_RewardedTrial.Choice(isnan(TimeStamp_RewardedTrial.Choice))=[];
TimeStamp_RewardedTrial.RewardRetrieval(isnan(TimeStamp_RewardedTrial.RewardRetrieval))=[];

% Unrewarded Trial
for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == -1
        Unreward_Trial(currentTrial) = currentTrial;
    else
        Unreward_Trial(currentTrial) = nan;
    end
    
end
Unreward_Trial(isnan(Unreward_Trial))=[];


for currentTrial= 1:SessionData.nTrials
    if SessionData.Reward(currentTrial) == -1
        TimeStamp_UnrewardedTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_UnrewardedTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_UnrewardedTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_UnrewardedTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_UnrewardedTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardedTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardedTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardedTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_UnrewardedTrial.CenterLitOn(isnan(TimeStamp_UnrewardedTrial.CenterLitOn))=[];
TimeStamp_UnrewardedTrial.Initiation(isnan(TimeStamp_UnrewardedTrial.Initiation))=[];
TimeStamp_UnrewardedTrial.Choice(isnan(TimeStamp_UnrewardedTrial.Choice))=[];
TimeStamp_UnrewardedTrial.RewardRetrieval(isnan(TimeStamp_UnrewardedTrial.RewardRetrieval))=[];

% Previous Reward history

for currentTrial= 2:SessionData.nTrials
    if SessionData.Reward(currentTrial-1) == 1
        TimeStamp_RewardAfterTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_RewardAfterTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_RewardAfterTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_RewardAfterTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_RewardAfterTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardAfterTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardAfterTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardAfterTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_RewardAfterTrial.CenterLitOn(isnan(TimeStamp_RewardAfterTrial.CenterLitOn))=[];
TimeStamp_RewardAfterTrial.Initiation(isnan(TimeStamp_RewardAfterTrial.Initiation))=[];
TimeStamp_RewardAfterTrial.Choice(isnan(TimeStamp_RewardAfterTrial.Choice))=[];
TimeStamp_RewardAfterTrial.RewardRetrieval(isnan(TimeStamp_RewardAfterTrial.RewardRetrieval))=[];

% Previous Unreward history

for currentTrial= 2:SessionData.nTrials
    if ~(SessionData.Reward(currentTrial-1) == 1)
        
        TimeStamp_UnrewardAfterTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_UnrewardAfterTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_UnrewardAfterTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_UnrewardAfterTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        
    else
        
        TimeStamp_UnrewardAfterTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardAfterTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardAfterTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardAfterTrial.RewardRetrieval(currentTrial) = nan;
        
    end
end

TimeStamp_UnrewardAfterTrial.CenterLitOn(isnan(TimeStamp_UnrewardAfterTrial.CenterLitOn))=[];
TimeStamp_UnrewardAfterTrial.Initiation(isnan(TimeStamp_UnrewardAfterTrial.Initiation))=[];
TimeStamp_UnrewardAfterTrial.Choice(isnan(TimeStamp_UnrewardAfterTrial.Choice))=[];
TimeStamp_UnrewardAfterTrial.RewardRetrieval(isnan(TimeStamp_UnrewardAfterTrial.RewardRetrieval))=[];

% Unrewarded(preRewarded)

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        ID_RewardAfterUnrewardTrial(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == 1 && ~(SessionData.Reward(currentTrial) == 1)
            ID_RewardAfterUnrewardTrial(currentTrial) = currentTrial;
        else
            ID_RewardAfterUnrewardTrial(currentTrial) = nan;
        end
    end
    
end
ID_RewardAfterUnrewardTrial(isnan(ID_RewardAfterUnrewardTrial))=[];

% Trial 1

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardAfterUnrewardTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == 1 && ~(SessionData.Reward(currentTrial) == 1)
            TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_RewardAfterUnrewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_RewardAfterUnrewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_RewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_RewardAfterUnrewardTrial.Choice(currentTrial) = nan;
            TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_RewardAfterUnrewardTrial.CenterLitOn(isnan(TimeStamp_RewardAfterUnrewardTrial.CenterLitOn))=[];
TimeStamp_RewardAfterUnrewardTrial.Initiation(isnan(TimeStamp_RewardAfterUnrewardTrial.Initiation))=[];
TimeStamp_RewardAfterUnrewardTrial.Choice(isnan(TimeStamp_RewardAfterUnrewardTrial.Choice))=[];
TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval(isnan(TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval))=[];

% Rewarded (previously unrewarded)
for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        ID_RewTrialPreUnrewarded(currentTrial) = nan;
    else
        if ~(SessionData.Reward(currentTrial-1) == 1) && SessionData.Reward(currentTrial) == 1
            ID_RewTrialPreUnrewarded(currentTrial) = currentTrial;
        else
            ID_RewTrialPreUnrewarded(currentTrial) = nan;
        end
    end
    
end
ID_RewTrialPreUnrewarded(isnan(ID_RewTrialPreUnrewarded))=[];


for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardAfterRewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardAfterRewardTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if ~(SessionData.Reward(currentTrial-1) == 1) && SessionData.Reward(currentTrial) == 1
            TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_UnrewardAfterRewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_UnrewardAfterRewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_UnrewardAfterRewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_UnrewardAfterRewardTrial.Choice(currentTrial) = nan;
            TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_UnrewardAfterRewardTrial.CenterLitOn(isnan(TimeStamp_UnrewardAfterRewardTrial.CenterLitOn))=[];
TimeStamp_UnrewardAfterRewardTrial.Initiation(isnan(TimeStamp_UnrewardAfterRewardTrial.Initiation))=[];
TimeStamp_UnrewardAfterRewardTrial.Choice(isnan(TimeStamp_UnrewardAfterRewardTrial.Choice))=[];
TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval(isnan(TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval))=[];

% Rewarded Trial (Previously Rewarded)(TrialLookup = 8)

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        ID_Reward_preReward(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == 1 && SessionData.Reward(currentTrial) == 1
            ID_Reward_preReward(currentTrial) = currentTrial;
        else
            ID_Reward_preReward(currentTrial) = nan;
        end
    end
    
end
ID_Reward_preReward(isnan(ID_Reward_preReward))=[];


for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_RewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_RewardAfterRewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_RewardAfterRewardTrial.Choice(currentTrial) = nan;
        TimeStamp_RewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == 1 && SessionData.Reward(currentTrial) == 1
            TimeStamp_RewardAfterRewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_RewardAfterRewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_RewardAfterRewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_RewardAfterRewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_RewardAfterRewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_RewardAfterRewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_RewardAfterRewardTrial.Choice(currentTrial) = nan;
            TimeStamp_RewardAfterRewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_RewardAfterRewardTrial.CenterLitOn(isnan(TimeStamp_RewardAfterRewardTrial.CenterLitOn))=[];
TimeStamp_RewardAfterRewardTrial.Initiation(isnan(TimeStamp_RewardAfterRewardTrial.Initiation))=[];
TimeStamp_RewardAfterRewardTrial.Choice(isnan(TimeStamp_RewardAfterRewardTrial.Choice))=[];
TimeStamp_RewardAfterRewardTrial.RewardRetrieval(isnan(TimeStamp_RewardAfterRewardTrial.RewardRetrieval))=[];

% Unrewarded Trial (Previously Unrewarded)(TrialLookup = 9)

for currentTrial= 1:SessionData.nTrials
    if currentTrial == 1
        TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_UnrewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
        TimeStamp_UnrewardAfterUnrewardTrial.Choice(currentTrial) = nan;
        TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
    else
        if SessionData.Reward(currentTrial-1) == -1 && SessionData.Reward(currentTrial) == -1
            TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
            TimeStamp_UnrewardAfterUnrewardTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
            TimeStamp_UnrewardAfterUnrewardTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
            TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
        else
            TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(currentTrial) = nan;
            TimeStamp_UnrewardAfterUnrewardTrial.Initiation(currentTrial) = nan;
            TimeStamp_UnrewardAfterUnrewardTrial.Choice(currentTrial) = nan;
            TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(currentTrial) = nan;
        end
    end
end

TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn(isnan(TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn))=[];
TimeStamp_UnrewardAfterUnrewardTrial.Initiation(isnan(TimeStamp_UnrewardAfterUnrewardTrial.Initiation))=[];
TimeStamp_UnrewardAfterUnrewardTrial.Choice(isnan(TimeStamp_UnrewardAfterUnrewardTrial.Choice))=[];
TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval(isnan(TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval))=[];

% BlockEnd Trial

for currentTrial= 1:SessionData.nTrials
    if SessionData.Tag.Trial_BlockEnd(currentTrial) == 1
        TimeStamp_BlockEndTrial.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_BlockEndTrial.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_BlockEndTrial.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_BlockEndTrial.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_BlockEndTrial.CenterLitOn(currentTrial) = nan;
        TimeStamp_BlockEndTrial.Initiation(currentTrial) = nan;
        TimeStamp_BlockEndTrial.Choice(currentTrial) = nan;
        TimeStamp_BlockEndTrial.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_BlockEndTrial.CenterLitOn(isnan(TimeStamp_BlockEndTrial.CenterLitOn))=[];
TimeStamp_BlockEndTrial.Initiation(isnan(TimeStamp_BlockEndTrial.Initiation))=[];
TimeStamp_BlockEndTrial.Choice(isnan(TimeStamp_BlockEndTrial.Choice))=[];
TimeStamp_BlockEndTrial.RewardRetrieval(isnan(TimeStamp_BlockEndTrial.RewardRetrieval))=[];

% BlockStart Trial
for currentTrial= 1:SessionData.nTrials
    if SessionData.Tag.Trial_BlockStart_1(currentTrial) == 1
        TimeStamp_BlockStartTrial_1.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_BlockStartTrial_1.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_BlockStartTrial_1.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_BlockStartTrial_1.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_BlockStartTrial_1.CenterLitOn(currentTrial) = nan;
        TimeStamp_BlockStartTrial_1.Initiation(currentTrial) = nan;
        TimeStamp_BlockStartTrial_1.Choice(currentTrial) = nan;
        TimeStamp_BlockStartTrial_1.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_BlockStartTrial_1.CenterLitOn(isnan(TimeStamp_BlockStartTrial_1.CenterLitOn))=[];
TimeStamp_BlockStartTrial_1.Initiation(isnan(TimeStamp_BlockStartTrial_1.Initiation))=[];
TimeStamp_BlockStartTrial_1.Choice(isnan(TimeStamp_BlockStartTrial_1.Choice))=[];
TimeStamp_BlockStartTrial_1.RewardRetrieval(isnan(TimeStamp_BlockStartTrial_1.RewardRetrieval))=[];
% BlockStart Trial_2
for currentTrial= 1:SessionData.nTrials
    if SessionData.Tag.Trial_BlockStart_2(currentTrial) == 1
        TimeStamp_BlockStartTrial_2.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_BlockStartTrial_2.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_BlockStartTrial_2.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_BlockStartTrial_2.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_BlockStartTrial_2.CenterLitOn(currentTrial) = nan;
        TimeStamp_BlockStartTrial_2.Initiation(currentTrial) = nan;
        TimeStamp_BlockStartTrial_2.Choice(currentTrial) = nan;
        TimeStamp_BlockStartTrial_2.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_BlockStartTrial_2.CenterLitOn(isnan(TimeStamp_BlockStartTrial_2.CenterLitOn))=[];
TimeStamp_BlockStartTrial_2.Initiation(isnan(TimeStamp_BlockStartTrial_2.Initiation))=[];
TimeStamp_BlockStartTrial_2.Choice(isnan(TimeStamp_BlockStartTrial_2.Choice))=[];
TimeStamp_BlockStartTrial_2.RewardRetrieval(isnan(TimeStamp_BlockStartTrial_2.RewardRetrieval))=[];



%% Trial Selection for PSTH
Trial_List = {'All' 'Rewarded','Unrewarded', ...
    'Trial_preRew', 'Trial_preUnrew', ...
    'UnrewTrial_preRew', 'RewaTrial_preUnrewd',...
    'RewTrial_preRew','UnrewaTrial_preUnrew',...
    'BlockEndTrial', 'BlockStartTrial_1', 'BlockStartTrial_2', 'Manual'
    };

% Manual input

for currentTrial= 1:SessionData.nTrials
    if sum(currentTrial == Manual_Trial) == 1
        TimeStamp_TEST.CenterLitOn(currentTrial) = TimeStamp_FP.CenterLitOn(currentTrial);
        TimeStamp_TEST.Initiation(currentTrial) = TimeStamp_FP.Initiation(currentTrial);
        TimeStamp_TEST.Choice(currentTrial) = TimeStamp_FP.Choice(currentTrial);
        TimeStamp_TEST.RewardRetrieval(currentTrial) = TimeStamp_FP.RewardRetrieval(currentTrial);
    else
        TimeStamp_TEST.CenterLitOn(currentTrial) = nan;
        TimeStamp_TEST.Initiation(currentTrial) = nan;
        TimeStamp_TEST.Choice(currentTrial) = nan;
        TimeStamp_TEST.RewardRetrieval(currentTrial) = nan;
    end
end

TimeStamp_TEST.CenterLitOn(isnan(TimeStamp_TEST.CenterLitOn))=[];
TimeStamp_TEST.Initiation(isnan(TimeStamp_TEST.Initiation))=[];
TimeStamp_TEST.Choice(isnan(TimeStamp_TEST.Choice))=[];
TimeStamp_TEST.RewardRetrieval(isnan(TimeStamp_TEST.RewardRetrieval))=[];


TrialLookup = 1;
% 
% % make  arrays based on behavior timestamps, then average rows into a mean vector for plotting

% [TrialLookup,tf] = listdlg('ListString',Trial_List);

if TrialLookup == 1
    TimeStamp_Choice_Analysis = TimeStamp_AllTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_AllTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_AllTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_AllTrial.RewardRetrieval;
elseif TrialLookup == 2
    TimeStamp_Choice_Analysis = TimeStamp_RewardedTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardedTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardedTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardedTrial.RewardRetrieval;
elseif TrialLookup == 3
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardedTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardedTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardedTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardedTrial.RewardRetrieval;
elseif TrialLookup == 4
    TimeStamp_Choice_Analysis = TimeStamp_RewardAfterTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardAfterTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardAfterTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardAfterTrial.RewardRetrieval;
elseif TrialLookup == 5
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardAfterTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardAfterTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardAfterTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardAfterTrial.RewardRetrieval;
elseif TrialLookup == 6
    TimeStamp_Choice_Analysis = TimeStamp_RewardAfterUnrewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardAfterUnrewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardAfterUnrewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardAfterUnrewardTrial.RewardRetrieval;
elseif TrialLookup == 7
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardAfterRewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardAfterRewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardAfterRewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardAfterRewardTrial.RewardRetrieval;
elseif TrialLookup == 8
    TimeStamp_Choice_Analysis = TimeStamp_RewardAfterRewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_RewardAfterRewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_RewardAfterRewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_RewardAfterRewardTrial.RewardRetrieval;
elseif TrialLookup ==9
    TimeStamp_Choice_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_UnrewardAfterUnrewardTrial.RewardRetrieval;
elseif TrialLookup == 10
    TimeStamp_Choice_Analysis = TimeStamp_BlockEndTrial.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_BlockEndTrial.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_BlockEndTrial.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_BlockEndTrial.RewardRetrieval;
elseif TrialLookup == 11
    TimeStamp_Choice_Analysis = TimeStamp_BlockStartTrial_1.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_BlockStartTrial_1.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_BlockStartTrial_1.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_BlockStartTrial_1.RewardRetrieval;
elseif TrialLookup == 12
    TimeStamp_Choice_Analysis = TimeStamp_BlockStartTrial_2.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_BlockStartTrial_2.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_BlockStartTrial_2.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_BlockStartTrial_2.RewardRetrieval;
elseif TrialLookup == 13
    TimeStamp_Choice_Analysis = TimeStamp_TEST.Choice;
    TimeStamp_Initiation_Analysis = TimeStamp_TEST.Initiation;
    TimeStamp_CenterLitOn_Analysis = TimeStamp_TEST.CenterLitOn;
    TimeStamp_RewardRetrieval_Analysis = TimeStamp_TEST.RewardRetrieval;
end

% Graph setting

% x-Axis Analysis window
nSecPrev = 10; %change to make different window
nSecPost = 10;

% x-Axis View window
nSecPrevWindow = 1; %change to make different window
nSecPostWindow = 1;

% y-axis
ymin = -1.5;
ymax = 1.5;

% convert seconds to TDT timestamps
nTsPrev = round (nSecPrev * samplingRate);
nTsPost = round (nSecPost * samplingRate);

totalTs = nTsPrev + nTsPost;
increment = (nSecPrev + nSecPost) / totalTs;
timeAxis = (-1 * nSecPrev) : increment : nSecPost;

% make raw_PSTH for Correct Press
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;


F(11) = figure('color', [1 1 1]); subplot(1,5,3)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Choice', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);

set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

legend ([h1 h2], 'pDMS', 'aDMS');

set(gca, 'layer', 'top');
set(gcf, 'Position', [0 600 1600 400])

%
% h1f=fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], 'g', 'EdgeColor', 'none');
% set(h1f,'facealpha',.5);
% h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color','r', 'Linewidth', 3);

% PsthArray_CORRECTPRESS_Block.Block1 = PsthArray_CORRECTPRESS(1:50, :);
% PsthArray_CORRECTPRESS_Block.Block2 = PsthArray_CORRECTPRESS(51:100, :);
% PsthArray_CORRECTPRESS_Block.Block3 = PsthArray_CORRECTPRESS(101:150, :);
% %PsthArray_CORRECTPRESS_Block.Block4 = PsthArray_CORRECTPRESS(151:200, :);



% make raw_PSTH for Initiation
nINITIATION = length(TimeStamp_Initiation_Analysis);
PsthArray_INITIATION_ch1 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_INITIATION_ch2 = NaN(nINITIATION,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nINITIATION
    thisTime = TimeStamp_Initiation_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_INITIATION_ch1(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_INITIATION_ch2(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_INITIATION.Err_INITIATION_ch1 = (nanstd(PsthArray_INITIATION_ch1))/sqrt(size(PsthArray_INITIATION_ch1,1));
FIGURE_INITIATION.Psth_INITIATION_ch1 = nanmean(PsthArray_INITIATION_ch1);
FIGURE_INITIATION.Err_Positive_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 + FIGURE_INITIATION.Err_INITIATION_ch1;
FIGURE_INITIATION.Err_Negative_ch1 = FIGURE_INITIATION.Psth_INITIATION_ch1 - FIGURE_INITIATION.Err_INITIATION_ch1;

FIGURE_INITIATION.Err_INITIATION_ch2 = (nanstd(PsthArray_INITIATION_ch2))/sqrt(size(PsthArray_INITIATION_ch2,1));
FIGURE_INITIATION.Psth_INITIATION_ch2 = nanmean(PsthArray_INITIATION_ch2);
FIGURE_INITIATION.Err_Positive_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 + FIGURE_INITIATION.Err_INITIATION_ch2;
FIGURE_INITIATION.Err_Negative_ch2 = FIGURE_INITIATION.Psth_INITIATION_ch2 - FIGURE_INITIATION.Err_INITIATION_ch2;




%F(12) = figure('color', [1 1 1]); 
subplot(1,5,2);
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch1, fliplr(FIGURE_INITIATION.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);


r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_INITIATION.Err_Positive_ch2, fliplr(FIGURE_INITIATION.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_INITIATION.Psth_INITIATION_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('Initiation', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 600 2000 400])


% PsthArray_INITIATION_Block.Block1 = PsthArray_INITIATION(1:50, :);
% PsthArray_INITIATION_Block.Block2 = PsthArray_INITIATION(51:100, :);
% PsthArray_INITIATION_Block.Block3 = PsthArray_INITIATION(101:150, :);
% PsthArray_INITIATION_Block.Block4 = PsthArray_INITIATION(151:200, :);

% make raw_PSTH for Liton
nCENTERLITON = length(TimeStamp_CenterLitOn_Analysis);
PsthArray_CENTERLITON_ch1 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CENTERLITON_ch2 = NaN(nCENTERLITON,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCENTERLITON
    thisTime = TimeStamp_CenterLitOn_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CENTERLITON_ch1(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CENTERLITON_ch2(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_CENTERLITON.Err_CENTERLITON_ch1 = (nanstd(PsthArray_CENTERLITON_ch1))/sqrt(size(PsthArray_CENTERLITON_ch1,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 = nanmean(PsthArray_CENTERLITON_ch1);
FIGURE_CENTERLITON.Err_Positive_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 + FIGURE_CENTERLITON.Err_CENTERLITON_ch1;
FIGURE_CENTERLITON.Err_Negative_ch1 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch1 - FIGURE_CENTERLITON.Err_CENTERLITON_ch1;

FIGURE_CENTERLITON.Err_CENTERLITON_ch2 = (nanstd(PsthArray_CENTERLITON_ch2))/sqrt(size(PsthArray_CENTERLITON_ch2,1));
FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 = nanmean(PsthArray_CENTERLITON_ch2);
FIGURE_CENTERLITON.Err_Positive_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 + FIGURE_CENTERLITON.Err_CENTERLITON_ch2;
FIGURE_CENTERLITON.Err_Negative_ch2 = FIGURE_CENTERLITON.Psth_CENTERLITON_ch2 - FIGURE_CENTERLITON.Err_CENTERLITON_ch2;



%F(13) = figure('color', [1 1 1]);
subplot(1,5,1)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch1, fliplr(FIGURE_CENTERLITON.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CENTERLITON.Err_Positive_ch2, fliplr(FIGURE_CENTERLITON.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CENTERLITON.Psth_CENTERLITON_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('CenterLightON', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 600 500 400])


% PsthArray_LITON_Block.Block1 = PsthArray_CENTERLITON(1:50, :);
% PsthArray_LITON_Block.Block2 = PsthArray_CENTERLITON(51:100, :);
% PsthArray_LITON_Block.Block3 = PsthArray_CENTERLITON(101:150, :);
% PsthArray_LITON_Block.Block4 = PsthArray_CENTERLITON(151:200, :);

nREWARD = length(TimeStamp_RewardRetrieval_Analysis);
PsthArray_REWARD_ch1 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_REWARD_ch2 = NaN(nREWARD,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nREWARD
    thisTime = TimeStamp_RewardRetrieval_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_REWARD_ch1(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_REWARD_ch2(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch2, thisIndex, nTsPrev, nTsPost);
end


FIGURE_REWARD.Err_REWARD_ch1 = (nanstd(PsthArray_REWARD_ch1))/sqrt(size(PsthArray_REWARD_ch1,1));
FIGURE_REWARD.Psth_REWARD_ch1 = nanmean(PsthArray_REWARD_ch1);
FIGURE_REWARD.Err_Positive_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 + FIGURE_REWARD.Err_REWARD_ch1;
FIGURE_REWARD.Err_Negative_ch1 = FIGURE_REWARD.Psth_REWARD_ch1 - FIGURE_REWARD.Err_REWARD_ch1;

FIGURE_REWARD.Err_REWARD_ch2 = (nanstd(PsthArray_REWARD_ch2))/sqrt(size(PsthArray_REWARD_ch2,1));
FIGURE_REWARD.Psth_REWARD_ch2 = nanmean(PsthArray_REWARD_ch2);
FIGURE_REWARD.Err_Positive_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 + FIGURE_REWARD.Err_REWARD_ch2;
FIGURE_REWARD.Err_Negative_ch2 = FIGURE_REWARD.Psth_REWARD_ch2 - FIGURE_REWARD.Err_REWARD_ch2;

%F(14) = figure('color', [1 1 1]);
subplot(1,5,4)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch1, fliplr(FIGURE_REWARD.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_REWARD.Err_Positive_ch2, fliplr(FIGURE_REWARD.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_REWARD.Psth_REWARD_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('RewardRetrieval', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([(-1 * nSecPrevWindow),nSecPostWindow]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [(-1 * nSecPrevWindow),nSecPostWindow], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [0 0], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
legend ([h1 h2], 'pDMS', 'aDMS');
legend('off')

set(gca, 'layer', 'top');
%set(gcf, 'Position', [500 0 500 400])
% PsthArray_REWARD_Block.Block1 = PsthArray_REWARD(1:50, :);
% PsthArray_REWARD_Block.Block2 = PsthArray_REWARD(51:100, :);
% PsthArray_REWARD_Block.Block3 = PsthArray_REWARD(101:150, :);
% PsthArray_REWARD_Block.Block4 = PsthArray_REWARD(151:200, :);

% make raw_PSTH for next ISI
nCHOICE = length(TimeStamp_Choice_Analysis);
PsthArray_CHOICE_ch1 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
PsthArray_CHOICE_ch2 = NaN(nCHOICE,nTsPrev+nTsPost+1); % preallocate arrays for speed
for i = 1:nCHOICE
    thisTime = TimeStamp_Choice_Analysis(i);
    thisIndex = round((thisTime*samplingRate))+1;
    PsthArray_CHOICE_ch1(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch1, thisIndex, nTsPrev, nTsPost);
    PsthArray_CHOICE_ch2(i,:) = processPhotDataRow_normDat(normDat_DBZ_ch2, thisIndex, nTsPrev, nTsPost);
end

FIGURE_CHOICE.Err_CHOICE_ch1 = (nanstd(PsthArray_CHOICE_ch1))/sqrt(size(PsthArray_CHOICE_ch1,1));
FIGURE_CHOICE.Psth_CHOICE_ch1 = nanmean(PsthArray_CHOICE_ch1);
FIGURE_CHOICE.Err_Positive_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 + FIGURE_CHOICE.Err_CHOICE_ch1;
FIGURE_CHOICE.Err_Negative_ch1 = FIGURE_CHOICE.Psth_CHOICE_ch1 - FIGURE_CHOICE.Err_CHOICE_ch1;

FIGURE_CHOICE.Err_CHOICE_ch2 = (nanstd(PsthArray_CHOICE_ch2))/sqrt(size(PsthArray_CHOICE_ch2,1));
FIGURE_CHOICE.Psth_CHOICE_ch2 = nanmean(PsthArray_CHOICE_ch2);
FIGURE_CHOICE.Err_Positive_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 + FIGURE_CHOICE.Err_CHOICE_ch2;
FIGURE_CHOICE.Err_Negative_ch2 = FIGURE_CHOICE.Psth_CHOICE_ch2 - FIGURE_CHOICE.Err_CHOICE_ch2;


%F(11) = figure('color', [1 1 1]); 
subplot(1,5,5)
hold on

o = 0.5;
r1 = 255;
g1 = 0;
b1 = 0;
rgb1_o = opacity (o, r1, g1, b1);
fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch1, fliplr(FIGURE_CHOICE.Err_Negative_ch1)], rgb1_o, 'EdgeColor', 'none');
h1= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch1,'Color',[(r1/255) (g1/255) (b1/255)], 'Linewidth', 3);

r2 = 0;
g2 = 0;
b2 = 255;
rgb2_o = opacity (o, r2, g2, b2);

fill([timeAxis, fliplr(timeAxis)],[FIGURE_CHOICE.Err_Positive_ch2, fliplr(FIGURE_CHOICE.Err_Negative_ch2)], rgb2_o, 'EdgeColor', 'none');
h2= plot (timeAxis,FIGURE_CHOICE.Psth_CHOICE_ch2,'Color',[(r2/255) (g2/255) (b2/255)], 'Linewidth', 3);

%labels, legend, make pretty, size
title ('ISI', 'FontSize', 10)
xlabel('Time (s)', 'FontSize', 20);
ylabel('z-score', 'FontSize', 20);

legend BOXOFF;
xlim ([0,8]);
xticks([0 3 8]);
ylim ([ymin, ymax]);
set(gca, ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'Fontsize', 18, ...
    'LineWidth'   , 2         );
set(gca, 'Color', 'none')
line('XData', [0,8], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')
line('XData', [3 3], 'YData', [ymin, ymax], 'LineStyle', '--', ...
    'LineWidth', 1, 'Color','k')

legend ([h1 h2], 'pDMS', 'aDMS');

set(gca, 'layer', 'top');
%set(gcf, 'Position', [0 0 2000 400])

