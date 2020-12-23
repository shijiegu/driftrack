
%% 
addpath(genpath('/home/lorenlab/Documents/MATLAB/Kilosort/')) % path to kilosort folder
addpath(genpath('/home/lorenlab/Documents/MATLAB/npy-matlab/')) % path to kilosort folder
addpath('/home/lorenlab/Documents/MATLAB/npy-matlab/') % for converting to Phy
rootZ = '/home/lorenlab/Documents/MATLAB/Kilosort/Shijie/NIckData/'; % the raw data binary file is in this folder
rootH = '/home/lorenlab/Documents/MATLAB/Kilosort/Shijie/NIckData/tmp_1ofsecond/';% path to temporary binary file (same size as data, should be on fast SSD)
pathToYourConfigFile = '/home/lorenlab/Documents/MATLAB/Kilosort/Shijie/NIckData/tmp_10ofsecond/'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'NP2_kilosortChanMap.mat';

%% A. load manipulator data
manip_positions = readNPY('manip.positions.npy');
manip_time_stamps = readNPY('manip.timestamps_p2.npy');
figure; plot(manip_time_stamps/60,manip_positions)

%% B. set time range
ops.trange    = [8*60 13*60]; % time range to sort
ops.NchanTOT  = 385; % total number of channels in your recording

run(fullfile(pathToYourConfigFile, 'configFile384.m'))
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);

%% C. Get Spike, Get Time, Location, Amplitude for all spikes.
fprintf('Looking for data inside %s \n', rootZ)

% main parameter changes from Kilosort2 to v2.5
ops.sig        = 20;  % spatial smoothness constant for registration
ops.fshigh     = 300; % high-pass more aggresively
ops.nblocks    = 5;   % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);
save(fullfile(rootH,'rez_step1.mat'),'rez','-v7.3')
%%
rez = datashift2(rez, 0); % last input is for shifting data
save(fullfile(rootH,'rez_uncorrected.mat'),'rez','-v7.3')

%% D. Plot spikes
del_t=10;
data=zeros(ceil(max(rez.st0(:, 2)))+1,ceil(max(rez.st0(:, 1)*del_t/ops.fs)));
for s=1:length(rez.st0(:, 2))
    l=ceil(rez.st0(s, 2))+1;
    t=ceil(rez.st0(s, 1)*del_t/ops.fs);
    data(l,t)=data(l,t)+rez.st0(s, 3);
end
% figure(194);
% imagesc(data);
% drawnow
% xlabel('time (0.1 s)')
% ylabel('spike position (um)')
% title('Drift map')
%filepath2save=['/home/lorenlab/Documents/drift_data/5_15min/data_',num2str(del_t),'of_a_second.mat'];
filepath2save=fullfile(rootH,'data.mat');
save(filepath2save,'data','-v7.3')
filepath2save(end-3:end)='.npy';
writeNPY(data,filepath2save);
%%
writeNPY(rez.dshift,fullfile(rootH,'benchmark_drift_10ofsecond.npy'));
writeNPY(rez.yblk,fullfile(rootH,'benchmark_y_10ofsecond.npy'));

%%
path='/home/lorenlab/Documents/MATLAB/Kilosort/Shijie/NIckData/tmp_10ofsecond/';
rez.dshift2=readNPY(fullfile(path,'pw4_drift_10ofsecond.npy'));
rez.dshift3=readNPY(fullfile(path,'r1pw4_drift_10ofsecond.npy'));

rez.ops.fproc3=rez.ops.fproc2;
rez.ops.fproc3(end-4)='3';

do_correction=1;
ops=rez.ops;
dshift=rez.dshift;
dshift2=rez.dshift2;
dshift3=rez.dshift3;
yblk=rez.yblk;
yblk2=rez.yblk2;
yblk3=rez.yblk2;
Nbatches=rez.temp.Nbatch;
if do_correction
    % sigma for the Gaussian process smoothing
    sig = rez.ops.sig;
    % register the data batch by batch
    dprev = gpuArray.zeros(ops.ntbuff,ops.Nchan, 'single');
    dprev2 = gpuArray.zeros(ops.ntbuff,ops.Nchan, 'single');
    dprev3 = gpuArray.zeros(ops.ntbuff,ops.Nchan, 'single');
    for ibatch = 1:rez.temp.Nbatch
        dprev = shift_batch_on_disk2(rez, ibatch, dshift(ibatch, :), yblk, sig, dprev, ops.fproc);
        dprev2 = shift_batch_on_disk2(rez, ibatch, dshift2(ibatch, :), yblk2, sig, dprev2, ops.fproc2);
        dprev3 = shift_batch_on_disk2(rez, ibatch, dshift3(ibatch, :), yblk3, sig, dprev3, ops.fproc3);
    end
    fprintf('time %2.2f, Shifted up/down %d batches. \n', toc, Nbatches)
else
    fprintf('time %2.2f, Skipped shifting %d batches. \n', toc, Nbatches)
end
%% 
figure(195);
st3=rez.st0;
set(gcf, 'Color', 'w')
% raster plot of all spikes at their original depths
st_shift = st3(:,2); %+ imin(batch_id)' * dd;
for j = 8:100
    % for each amplitude bin, plot all the spikes of that size in the
    % same shade of gray
    ix = st3(:, 3)==j; % the amplitudes are rounded to integers
    plot(st3(ix, 1)/ops.fs, st_shift(ix), '.', 'color', [1 1 1] * max(0, 1-j/40)) % the marker color here has been carefully tuned
    hold on
end
axis off
AxesH = gca;
drawnow;
InSet = get(AxesH, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
%%
Outer_pos=get(gca, 'OuterPosition');
disp(Outer_pos)
set(AxesH, 'Position', [0, 0, Outer_pos(3),Outer_pos(4)])

%%
f=getframe(gcf);
set(gcf,'units','points','position',[50,50,300,400])
data=mean(f.cdata,3);

%xlabel('time (sec)')
%ylabel('spike position (um)')
%title('Drift map')


%% ONLY DO THIS AFTER YOU HAVE RUN THE OLD METHOD, HACKISH WAY HERE
old_whitened=rez.ops.fproc;
new_whitened=rez.ops.fproc2;
new_whitened_r1=rez.ops.fproc3;
%%
rez.ops.fproc=new_whitened_r1;
rez.ops.fproc2=old_whitened;
%%
rez.ops.fproc='/home/lorenlab/Documents/MATLAB/Kilosort/Shijie/NIckData/tmp_10ofsecond/temp_wh_backup.mat';
%%
rootZ2=fullfile(rootH,'unwarped');
mkdir(rootZ2)
% ORDER OF BATCHES IS NOW RANDOM, controlled by random number generator
iseed = 1;
                 
% main tracking and template matching algorithm
rez = learnAndSolve8b(rez, iseed);

% OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
% See issue 29: https://github.com/MouseLand/Kilosort/issues/29
rez = remove_ks2_duplicate_spikes(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% decide on cutoff
rez = set_cutoff(rez);
% eliminate widely spread waveforms (likely noise)
rez.good = get_good_units(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ2);

%% if you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

% Ensure all GPU arrays are transferred to CPU side before saving to .mat
rez_fields = fieldnames(rez);
for i = 1:numel(rez_fields)
    field_name = rez_fields{i};
    if(isa(rez.(field_name), 'gpuArray'))
        rez.(field_name) = gather(rez.(field_name));
    end
end

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(rootZ, 'rez2.mat');
save(fname, 'rez', '-v7.3');