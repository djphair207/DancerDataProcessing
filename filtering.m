close all

filename = "AbbyW_01.mat";
TPS = 100;

data = importdata(filename);                                  % read in the accelerometer data from file
t = transpose(0:data.x_values.increment:data.x_values.increment*(data.x_values.number_of_values-1));% this only works for the data collect from the SIEMENS system
numTaps = fix(t(end)*TPS);
A0 = data.y_values.values(:,1);
A1 = data.y_values.values(:,2);
SampFreq = size(data.y_values.values,1)/t(end);

%BPfilt = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',10,'HalfPowerFrequency2',2000,'SampleRate',SampFreq);     % create a bandpass filter
%BPfilt = designfilt('bandpassiir','FilterOrder',8,'HalfPowerFrequency1',10,'HalfPowerFrequency2',2000,'SampleRate',SampFreq);     % create a bandpass filter

BPfilt = designfilt('bandpassiir','FilterOrder',6,'HalfPowerFrequency1',200,'HalfPowerFrequency2',1000,'SampleRate',SampFreq);     % create a bandpass filter


%BPfilt = designfilt('bandpassiir','FilterOrder',4,'CutoffFrequency1',10,'CutoffFrequency2',2000,'SampleRate',SampFreq);     % create a bandpass filter
%BPfilt = designfilt('bandpassfir', 'FilterOrder', 4, 'CutoffFrequency1', 100, 'CutoffFrequency2', 500, 'SampleRate',SampFreq); 

filt_A0 = filtfilt(BPfilt,A0);     % filter both arrays of accelerometer data 
filt_A1 = filtfilt(BPfilt,A1);

%% Fourier Transform
% 
% fft_A0 = fft(A0);
% fft_A1 = fft(A1);
% 
% 

%% plotting
figure;
hold on
plot(t,filt_A0,"r")
plot(t,filt_A1,"b")
title('Accelerometer Response');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
hold off

% figure;
% hold on
% plot(t,A0,"r")
% plot(t,A1,"b")
% title('Accelerometer Response');
% xlabel('Time (s)');
% ylabel('Acceleration (m/s^2)');
% hold off