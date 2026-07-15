%% Homework Imaging for Neuroscience - second part of the course
% Homework 4 - second part (last 3 CFU)
% Group components: Lia Frizzarin, Michele Simoncelli, Tommaso Voltolina

% The aim is to compare DOT reconstructed images obtained using different 
% lambdas for the regularization when solving the inverse problem. 
% DOT data acquired in one adult human performing a color-naming task.

clearvars
close all
clc

%% preparing the paths
basePath = pwd;

addpath(genpath(fullfile(basePath,"iso2mesh")))
addpath(genpath(fullfile(basePath,"homer2")))

%% loading NIRS data
load('S126_Color.nirs','-mat');

nCh = size(SD.MeasList,1)/2; % number of channels

fs = 1/(t(2)-t(1)); % compute sampling frequency


%% 1. Plot the 3D array configuration
figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',12) %3D plot of sources
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.g','MarkerSize',12) %3D plot of detectors
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'k') %3D plot of channels
end
legend('Sources', 'Detectors', 'Channels',Location='northeast')

% overlapping channels 

%% 2. Source-detector distance 
% Compute Source-Detector (SD) distance for each channel
distCh = zeros(nCh,1);
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    distCh(iCh) = sqrt(sum((src-det).^2));
end

% Plot the SD distance as histogram
figure;
histogram(distCh,20)
xlabel('SD distance [mm]')
ylabel('N of channels')

% channels of different lengths


%% 3. Identify "bad" channels
% Remove noisy channels
dRange = [0.03 3];
SNRrange = 7;
remCh = removeNoisyChannels(d,dRange,SNRrange);

% place the remCh vector in the SD struct so that the functions 
% will know which channel is good and should be processed and which not
SD.MeasListAct = remCh; 

% Plot the 3D array configuration highlighting the bad channels with a different color 
figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',12) %3D plot of sources
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.g','MarkerSize',12) %3D plot of detectors
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    if remCh(iCh) == 0
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'c') %3D plot of bad channels
    else
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'k') %3D plot of good channels
    end
end
xlabel('x')
ylabel('y')
zlabel('z')
legend('Sources', 'Detectors', 'Channels',Location='northeast')

% Plot the array configuration keeping in the plot only the good channels. 
GoodChannels = find(SD.MeasListAct == 1);
figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',12) %3D plot of sources
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.g','MarkerSize',12) %3D plot of detectors
for iCh = 1:length(GoodChannels)
    src = SD.SrcPos(SD.MeasList(GoodChannels(iCh),1),:);
    det = SD.DetPos(SD.MeasList(GoodChannels(iCh),2),:);
    plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'k') %3D plot of channels
end
legend('Sources', 'Detectors', 'Good Channels',Location='northeast')


%% 4. Pre-processing steps of fNIRS data
%% 4.a Convert to optical density changes
meanValue = mean(d);
dodConv = -log(abs(d)./meanValue);

% Visualize optical density changes of non-removed channels
dodDataRemoved = dodConv;
dodDataRemoved(:,remCh==0) = [];

figure; 
subplot(211)
plot(t,dodDataRemoved(:,1:end/2))
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')
subplot(212)
plot(t,dodDataRemoved(:,end/2+1:end))
title('Second Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

%% 4.b Wavelet motion artifacts detection and correction
% Set iqr to 0.5
iqr = 0.5;
% Run wavelet motion correction...
%dodWavelet = hmrMotionCorrectWavelet(dodConv,SD,iqr);
% ...or load the already computed dodWavelet.mat
load("dodWavelet.mat")


%% 4.c Band-pass filter
lowerCutOff = 0.01;
higherCutOff = 0.5;
dodFilt = hmrBandpassFilt(dodWavelet,fs,lowerCutOff,higherCutOff);

% checking the goodness of the band-pass filter
dodWavGood = dodWavelet(:,remCh==1);
dodFiltGood = dodFilt(:,remCh==1);
figure;
plot(t,dodWavGood(:,1))
hold on;
plot(t,dodFiltGood(:,1))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')
legend('original','band-pass filtered')


%% 4.d Compute block-averaged hemodynamic response on the optical density data
tRange = [-2 15]; % range of timimg around stimulus to define a trial
sRange = fix(tRange*fs); % convert the time in seconds to samples

tHRF = tRange(1):1/fs:tRange(2); % time vector for the hemodynamic response (and trials)

% initialize the matrix that will contain our average hemodynamic response for each channel (for both wavelength) and condition
dodAvg = zeros(length(tHRF),size(dodFilt,2),size(s,2));

for iS = 1:size(s,2) % for each condition

    % Get the timing of stimulus presentation for that condition
    stimulusTiming = find(s(:,iS)==1); 

    % Initialize the matrix that will contain the single trial responses
    % for that condition for all channels at both wavelengths
    ytrial = zeros(length(tHRF),size(dodFilt,2),length(stimulusTiming));
    
    nTrial = 0;
    for iT = 1:length(stimulusTiming) % for each stimulus presented (for eacht trial)
        if (stimulusTiming(iT)+sRange(1))>=1 && (stimulusTiming(iT)+sRange(2))<=size(dodFilt,1) % Check that there are enough data pre and post stimulus (this is useful to check that the first stimulus is presented at least 2 seconds after the start of the acquisition and that the last stimulus has at least 18 seconds of data afterwards)
            nTrial = nTrial + 1;
            ytrial(:,:,nTrial) = dodFilt(stimulusTiming(iT)+[sRange(1):sRange(2)],:); % extract the trial from the dc data
        end
    end
    
    % Average trials (the fourth dimension of the ytrial matrix)
    dodAvg(:,:,iS) = mean(ytrial(:,:,1:nTrial),3);

    % Correct for the baseline
    for ii = 1:size(dodAvg,2) % for each channel
        foom = mean(dodAvg(1:-sRange(1),ii,iS),1); % compute baseline as average of the signal in the -2:0 seconds time range
        dodAvg(:,ii,iS) = dodAvg(:,ii,iS) - foom; % subtract the baseline from the average hemodynamic responses
    end
end


%% 5. Head volume mesh 

% Load data
load(fullfile('MNI152_headModel','HeadVolumeMesh.mat'))

% Plot volumetric head mesh and, overlaid on that, the sources and detectors of the array
figure;
hold on
plotmesh(HeadVolumeMesh.node(:,1:3),HeadVolumeMesh.elem)
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',15)
axis("vis3d")
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.g','MarkerSize',15)
axis("vis3d")
legend('Mesh', 'Sources', 'Detectors',Location='northeast')


%% 6. GM Volumetric Mesh
% Display the whole array sensitivity for the first wavelength on the volumetric GM mesh with all channels
% and, in a separate figure, by removing the “bad” channels

% Load Jacobian
load('GRINP.jac','-mat')

HeadVolumeMesh.node(:,4) = (sum(J{1}.vol));

% Display whole array sensitivity on GM volume mesh for all the channels
figure;
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
caxis([-3 0])
title("whole array sensitivity, first wavelength, all channels")

% Remove bad channels from Jacobian
JGood = J;
for i = 1:length(SD.Lambda) % For each wavelength (we have two Js)
    tmp = J{i}.vol;
    JGood{i}.vol = tmp(SD.MeasListAct(SD.MeasList(:,4)==i)==1,:);
end
HeadVolumeMesh.nodeGood = HeadVolumeMesh.node;
HeadVolumeMesh.nodeGood(:,4) = (sum(JGood{1}.vol));

% Display whole array sensitivity on GM volume mesh for the 'good' channels
% only
figure;
plotmesh(HeadVolumeMesh.nodeGood,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
caxis([-3 0])
title("whole array sensitivity, first wavelength, only good channels")


%% 7. Reconstruct HbO and HbR images mapped to the surface GM mesh

% Jacobian containing only the good channels
JCropped = {JGood{1}.vol, JGood{2}.vol};

% we want to compare the results of the inverse problem solved for
% different values of the regularization parameter
regularization = 1; % = [1,2,3,4,5]

switch regularization
    case 1
        lambda1 = 0.0001;  % almost no regularization
    case 2
        lambda1 = 0.01;    % weak regularization
    case 3
        lambda1 = 0.1;     % normal regularization
    case 4
        lambda1 = 1;       % strong regularization
    case 5
        lambda1 = 10;      % over-regularization
end

% Compute inverse of Jacobian
invJ = cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) %for each Jacobian
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    S=svd(JJT);
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda1*max(S)));
end

% Data to reconstruct are optical density changes compared to a baseline.
% In our case the baseline is 0, therefore we want to reconstruct 0-our
% data
datarecon = -dodAvg;

% Inizialize matrices and load useful stuff
load(fullfile('MNI152_headModel','GMSurfaceMesh.mat')) % Load GM surface and surface mesh
nNodeVol = size(HeadVolumeMesh.node,1);  % The node count of the volume mesh
nNodeGM = size(GMSurfaceMesh.node,1); % The node count of the GM mesh
nFrames = size(datarecon,1); % Number of samples to reconstruct
load('vol2gm.mat')
wavelengths = SD.Lambda; % wavelengths of the system
nWavs = length(SD.Lambda); % n of wavelengths

% Initialize final results matrices
hbo.vol = zeros(nFrames,nNodeVol);
hbr.vol = zeros(nFrames,nNodeVol);
hbo.gm = zeros(nFrames,nNodeGM);
hbr.gm = zeros(nFrames,nNodeGM);

% Obtain specific absorption coefficients
Eall = [];
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    Etmp = Etmp(1:2); % HbO and HbR only
    Eall = [Eall; Etmp./1e7]; % This will be nWavs x 2;
end

% Perform reconstruction for each frame with the multispectral image
% reconstruction approach
    
% For each frame
for frame = 1:nFrames

    % Reconstruct absorption changes
    muaImageAll = zeros(nWavs,nNodeVol);
    for wav = 1:nWavs
        dataTmp = squeeze(datarecon(frame,SD.MeasList(:,4)==wav & SD.MeasListAct==1));
        invJtmp = invJ{wav};
        tmp = invJtmp * dataTmp';
        muaImageAll(wav,:) = tmp; % This will be nWavs * nNode
    end

    % or
    Eallinv = pinv(Eall);
    tmp = Eallinv*muaImageAll;
    hbo_tmpVol = tmp(1,:);
    hbr_tmpVol = tmp(2,:);

    % Map to GM surface mesh
    hbo_tmpGM = (vol2gm*hbo_tmpVol');
    hbr_tmpGM = (vol2gm*hbr_tmpVol');

    % Book-keeping and saving
    hbo.vol(frame,:) = hbo_tmpVol;
    hbr.vol(frame,:) = hbr_tmpVol;
    hbo.gm(frame,:) = hbo_tmpGM;
    hbr.gm(frame,:) = hbr_tmpGM;

end

% Plot reconstructed images
tRecon = [0; 7; 15]; % time points 0s, 7s, 15s
baseline = abs(tRange(1)); % two seconds of baseline

for ii = 1:length(tRecon) 
    % select the time point
    tt = tRecon(ii);

    sRecon = fix(tt*fs)+fix(baseline*fs); % Convert to samples
    load greyJet % load colormap to make better images
    
    % Assign image to fourth column of node
    GMSurfaceMesh.node(:,4) = hbo.gm(sRecon,:);
    figure;
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    caxis([-0.05 0.05]) % Set the limit of the colorbar
    view([0 90]) % Set the view angle
    title(['HbO, t=' num2str(tt) 's'])
    colormap(greyJet) % set the loaded colormap
    hb = colorbar;
    hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
    axis off % remove axis
    
    GMSurfaceMesh.node(:,4) = hbr.gm(sRecon,:);
    figure;
    plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
    view([0 90])
    caxis([-0.05 0.05])
    title(['HbR, t=' num2str(tt) 's'])
    colormap(greyJet)
    hb = colorbar;
    hb.Label.String = {'\DeltaHbR [\muM]'};
    axis off
end


%% 8. Comparing the results obtained with different regularization parameters

% lambda1=0.0001 --> we are detecting strong spots of activation for both
%                    HbR and HbO, at all the three time points; we are almost surely
%                    overfitting the data, thus the solution contains also
%                    the noise.
% lambda1=0.01   --> we are able to detect quite strong activation spots for
%                    both HbR and HbO, however, we run an higher risk of
%                    overfitting the data.
% lambda1=0.1    --> it is the standard values of the regularization
%                    parameter, we are able to detect activation spots for both HbR (weaker)
%                    and HbO (stronger), eapecially for time point 7s.
% lambda1=1      --> we are oversmoothing the signal; we are able to detect
%                    only weak activation spots for HbO, which tipically reaches higher values
%                    than HbR but it is more contaminated by physiological
%                    noise.
% lambda1=10     --> we are oversmoothing the signal; we cannot see any
%                    activation spot.
%
% The optimal lambda1 may lay in the 0.01-0.1 range.


