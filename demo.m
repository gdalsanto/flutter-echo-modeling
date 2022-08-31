clear variables; close all; clc

% enclosing directory
mainFolder = pwd; 

% include auxiliary directories
dirList = ["functions", "audio", "data", fullfile("audio","RIRs"), fullfile("audio","output")];
for i = 1:length(dirList)
    addpath(fullfile(pwd, dirList(i)))
end

%% 1 Reference Room Impulse Response

% load rir containing flutter 
[refRIR, sampleRate] = audioread("rir_arni_c16c8m6s1.wav");
refRIR = refRIR/max(abs(refRIR));

% load diffuse rir
diffRIR = audioread("rir_arni_c16c12m6s1.wav");
diffRIR = diffRIR/max(abs(diffRIR));

% allign signals 
sampleDelay = finddelay(refRIR,diffRIR);
if sampleDelay
    if sampleDelay >= 0
        diffRIR = diffRIR(sampleDelay:end);
        nPad = length(refRIR)-length(diffRIR);
        diffRIR = [diffRIR; zeros(nPad,1)];
    else
        refRIR = refRIR(abs(sampleDelay):end);
        nPad = length(diffRIR)-length(refRIR);
        refRIR = [refRIR; zeros(nPad,1)];
    end
end

% ----- 1.1 Find onset sample ----- %
winSize = 441;  % initial segment length 
initEnergy = 10*log10(refRIR(1:winSize).*refRIR(1:winSize));

% set energy theshold 
thVal = mean(initEnergy(~isinf(initEnergy))) + 20*log10(20);  

% find onset sample
initSilence = find(10*log10(refRIR.*refRIR)>thVal,1,'first');

% ----- 1.2 Expected repetition time ----- %
c = 343;        % speed of sound
lenRoom = 8.1;  % length of variable acoustics room (Arni) in meters
expRepTonality = c/lenRoom/2;       % expected repetition tonality from room geometry 
expRepTime = 1/expRepTonality;      % expected repetition time t_r

%% 2 Flutter echo parameters

% ----- 2.1 repetition time ----- %
[repetitionTime, flutteryRange, lateRevTime] = repetition_time('rir_arni_c16c8m6s1.wav', 4, expRepTime, 0);
repetitionTime = repetitionTime*10;
% derive decay slope simulated by the attenuation filter
% load the filter coefficents from the data folder
varNames = ["aH.mat", "aL.mat", "bL.mat", "bH.mat", "peakRT.mat", "kMax.mat"];
for i = 1 : length(varNames)
    load(fullfile(dirList(3), varNames(i)));
end
peakRT = peakRT-0.105;

% filter design parameters 
repTimeSamples = round(sampleRate*repetitionTime);
L = repTimeSamples;                     % delay length
R = round(1/kMax*L);                    % ripple delay length
K = 5;                                  % downsampling factor 

refGain = -15;                          % reference gain in dB
nfft = round(repetitionTime*sampleRate)*2;  % number of fft points

hGeqH = impz(bH, aH, nfft);
hGeq = filter(bL, aL, hGeqH);
hGeq = hGeq*10^(peakRT/20);

% filter with ripple filter
r = 0.02;
hLoop = filter([r zeros(1, R-1) -1], 1, hGeq);
[HLoop, ~] = freqz(hLoop, 1, nfft, sampleRate);
maxHLoop = max(mag2db(abs(HLoop))); 
% match the peak gain
hLoop = hLoop.*10^((peakRT-maxHLoop)/20);
[HLoop, ~] = freqz(hLoop, 1, nfft, sampleRate);
loopRT = db2rt(db(HLoop),L,sampleRate);
slopeRT = -60./loopRT;
slopeRT = [slopeRT(1:end-1); slopeRT(end-1:-1:1)];

% ----- 2.2 generate intial pulse ----- %
initPulse = pulse_shape(refRIR, repetitionTime, lateRevTime, initSilence, repetitionTime, slopeRT, sampleRate, 0);
% resample initial pulse shape
initPulse = resample(initPulse,1,K);
L = length(initPulse); 
R = round(1/kMax*L); 

density = round(2205/K);            % impulses per second 
Td = round(sampleRate/K/density);   % grid size (length of EVN) 
delta = 1;                          % range where the impulses occur in samples

% generate velvet noise sequence 
[~, vns, ~] = velvetnsmkx(L, density, delta, sampleRate/K);

% take the envelope of the inital pulse shape
[envUp, envLow] = envelope(initPulse/max(abs(initPulse)),5,'peak');
% normalize the envelope 
envUp = envUp/max(abs(envUp));
envLow = envLow/max(abs(envLow));

buffVns = vns;
vns(buffVns>0) = buffVns(buffVns>0).*envUp(buffVns>0)';
vns(buffVns<0) = -buffVns(buffVns<0).*envLow(buffVns<0)';


%% 3 Flutter echo synthesis

% ----- 3.1 Read input signal ----- %
[inputSig, ~] = audioread("speech.wav");
inputSig = inputSig(:,1)/max(abs(inputSig(:,1)));    
inputSig = inputSig(1:end-mod(length(inputSig),K));

% add trailing zeros for reverb tail
nSamples = 4*sampleRate + max(size(inputSig));
inputSig = [inputSig(:)', zeros(1, nSamples-length(inputSig))];

% unit impulse of sample length
unitImpulse = [zeros(1, initSilence), 1, zeros(1, nSamples-initSilence-1)];

% ----- 3.2 Downsampling ----- % 
% antialiasing filter
beta = 5; 
order = 10; 
decFilter = dsp.FIRDecimator(K, fir1(order,1/K,kaiser(order+1,beta)));
inputSigK = decFilter(inputSig');

% adjust parameters according to new sample rate
initSilenceK = round(initSilence/K);
nSamplesK = round(nSamples/K);
unitImpulseK = [zeros(1, initSilenceK), 1, zeros(1, nSamplesK-initSilenceK-1)];


% ----- 3.3 Attenuation filter design ----- %
% load the filter coefficents from the data folder
varNames = ["aHdown.mat", "aLdown.mat", "bLdown.mat", "bHdown.mat"];
for i = 1 : length(varNames)
    load(fullfile(dirList(3), varNames(i)));
end
nfft = 8192; 

hGeqH = impz(bHdown, aHdown, nSamplesK/4);
hGeq = filter(bLdown, aLdown, hGeqH);
hGeq = hGeq*10^(peakRT/20);

% filter with ripple filter
hLoop = filter([r zeros(1, R-1) -1], 1, hGeq);
[HLoop, freqAxisPlot] = freqz(hLoop, 1, nfft, sampleRate/K);
maxHLoop = max(mag2db(abs(HLoop))); 
% match the peak gain
hLoop = hLoop.*10^((peakRT-maxHLoop)/20);
% extract group delay 
[groupDelay, ~] =  grpdelay(hLoop,1,nfft); 
% it can happen that grpdelay gives high results (>1e6)
avrGroupDelay = mean(groupDelay(abs(groupDelay)<1e6)); 
% compensate for group delay 
L = L - round(avrGroupDelay);
% filters in feedforward branch
cir = filter([zeros(1, (L-2)) hLoop'],[1 zeros(1,(L-2)) hLoop'], inputSigK);

% ----- 3.4 Velvet noise sequence convolution ----- %
% velvet noise convolution
lateRevK = filter(vns, 1, cir);

% the following is just to then adjust the power more easily 
flutterFiltK = filter([zeros(1, (L-2)) hLoop'],[1 zeros(1,(L-2)) hLoop'], unitImpulseK);
flutterFiltK = filter(vns, 1, flutterFiltK);

% ----- 3.5 Upsampling ----- % 
order = 20;
beta = 5;
intFilter = dsp.FIRInterpolator(K, fir1(order,1/K,kaiser(order+1,beta)));
lateRev = intFilter(lateRevK)';
flutterFilt = intFilter(flutterFiltK')';

% ----- 3.6 Power Adjustment ----- % 
% reference RIR 
sampleRange = round(repTimeSamples*1.2):length(refRIR);
refPower = compute_power(refRIR, sampleRange, 0);
% diffuse RIR 
earlyPower = compute_power(diffRIR, sampleRange, 0);
% synthesized RIR
flutterPower = compute_power(flutterFilt, sampleRange, 0);

% compute gain for power adjustment
targetPower = (sqrt(refPower) - sqrt(earlyPower))^2;
gain = sqrt(targetPower/flutterPower);

lateRev = lateRev*gain;

%% 4 Result
% generate early reflections
earlyRef = filter(diffRIR(initSilence:end),1,inputSig);

wetGain = 15;
synthRIR = earlyRef  + lateRev*wetGain; 
audiowrite(fullfile(dirList(end),'output_10tr.wav'),synthRIR/max(abs(synthRIR)),sampleRate);
