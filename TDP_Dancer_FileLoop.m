function result = TDP_Dancer_FileLoop(varargin)
% TDP_Dancer_FileLoop : Reads data from tensiometer data files, filters,
% interpolates, calculates wavespeed at each point, plots wavespeed
% vs time

close all
%disp(1)
%% Set up for varargin
% For a more complete description of how this works, see the MATLAB
% documentation for 'varargin'

defaultTPS = 100;                % default taps per second is 50 - THIS WAS CHANGED FOR THE DATA FROM THE DANCERS COLLECTED IN NOV 2023
defaultMinSpeed = 5;            % default Minimum wave speed is 5 m/s
defaultFilterOrder = 4;         % default filter order is 4
defaultCutFreq = 1;             % default cut off frequency is 1 Hz
defaultInterpAmount = 100;      % default number of interpolation points is 100
defaultInterpMethod = 1;        % 1 - interp1, 2 - fft - DON'T USE FFT YET
defaultSpacing = 10;            % default distance between accelerometers is 10 mm

p = inputParser;                                                                    % create parser object
validTPS = @(x) isnumeric(x) && (x > 0);                                            % define valid inputs
validMinSpeed = @(x) isnumeric(x) && (x > 0);
validFilterOrder = @(x) isnumeric(x) && (0 < x) && (x < 20) && isscalar(x);
validCutFreq = @(x) isnumeric(x) && (x > 0);
validInterpAmount = @(x) isnumeric(x) && (x > 0);
validInterpMethod = @(x) (x == 1) || (x==2);
validSpacing = @(x) isnumeric(x) && (x > 0);

% Add objects to the parser
addParameter(p,"TPS",defaultTPS,validTPS);
addParameter(p,"MinSpeed",defaultMinSpeed,validMinSpeed);
addParameter(p,"FilterOrder",defaultFilterOrder,validFilterOrder);
addParameter(p,"CutoffFreq",defaultCutFreq,validCutFreq);
addParameter(p,"InterpAmount",defaultInterpAmount,validInterpAmount);
addParameter(p,"InterpMethod",defaultInterpMethod,validInterpMethod);
addParameter(p,"Spacing",defaultSpacing,validSpacing);

parse(p,varargin{:});                                                      % Parse the user input
%disp(2)

%% accelerometer data
n = 1;
while n <= 25
    filename = uigetfile;
    data = importdata(filename);                                  % read in the accelerometer data from file
    t = transpose(0:data.x_values.increment:data.x_values.increment*(data.x_values.number_of_values-1));% this only works for the data collect from the SIEMENS system
    numTaps = fix(t(end)*p.Results.TPS);
    A0 = data.y_values.values(:,1);
    A1 = data.y_values.values(:,2);
    SampFreq = size(data.y_values.values,1)/t(end);
    interp_t = linspace(t(1), t(end), size(t,1)*p.Results.InterpAmount);    % make a new time vector based on interp amount 
    TimeToIndexConversion = length(interp_t)/max(t);                        % this is the value that converts from 'time space' to 'index space' & vice versa
    
    %disp(3)
    
    %% Filtering
    BPfilt = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',200,'HalfPowerFrequency2',1500,'SampleRate',SampFreq);     % create a bandpass filter
    filtered_original_A0 = filtfilt(BPfilt,A0);
    filtered_original_A1 = filtfilt(BPfilt,A1);
    
    %% Interpolate the new y and t values using interp1
    
    filtered_A0 = interp1(t,filtered_original_A0,interp_t);
    filtered_A1 = interp1(t,filtered_original_A1,interp_t);
    
    %% Calculating where each tap from the tapper will occur
    
    % / / / / / testing new triggering method / / / / / %
    
    % instead of just starting at the beginning of the line and looking for a
    % tap, let's try windowing a section where we KNOW there will be a tap. To
    % do this, divide the length of the interpolated accel data vector by the
    % number of taps. This should guarantee that a tap will be found in the
    % window.
    
    windowMax = length(interp_t)/numTaps;
    
    x = zeros(1,numTaps);
    i = 1;
    j = 1;
    foundTap = false;
    
    window = filtered_A0(1:windowMax);
    windowAvg = sum(window)/length(window);
    windowStd = std(window);
    
    while foundTap == false
        i = i+1;
        if filtered_A0(i) > windowAvg + 3*windowStd || filtered_A0(i) < windowAvg - 3*windowStd
            foundTap = true;
        end
    end
    
    if foundTap == true
        while i < length(interp_t)
            x(j) = fix(i);
            j = j+1;
            i = fix(i + (1/p.Results.TPS)*TimeToIndexConversion);
        end
    else
        disp("foundTap was not marked true")
    end
    % / / / / / / / / / / / / / %
    
    %% calculating the speed of the signal at each tap
    MaxLags = fix(.01*TimeToIndexConversion/p.Results.MinSpeed); % MinSpeed -> PARAMETER
    Lower_bound = q*10*p.Results.InterpAmount;
    Upper_bound = (60+(m*10))*p.Results.InterpAmount;
    filtered_interp1_speed = zeros(1,numTaps);
    time = zeros(1,numTaps);
    filtered_difference = zeros(1,numTaps);
    for i = 1:length(x)%taps      % walk through all the taps, windowing each of them
        if x(i)-Lower_bound < 0     % if the window extends to neg time, use 0 instead
            windowed_time = interp_t(1:x(i) + Upper_bound);
            new_filtered_A0 = filtered_A0(1:x(i) + Upper_bound);
            new_filtered_A1 = filtered_A1(1:x(i) + Upper_bound);
        elseif x(i) + Upper_bound > x(end)  % if the window extend outside the big end of x, use the last x val instead
            windowed_time = interp_t(x(i) - Lower_bound:x(end));
            new_filtered_A0 = filtered_A0(x(i) - Lower_bound:x(end));
            new_filtered_A1 = filtered_A1(x(i) - Lower_bound:x(end));
        else
            windowed_time = interp_t(x(i) - Lower_bound:x(i) + Upper_bound);
            new_filtered_A0 = filtered_A0(x(i) - Lower_bound:x(i) + Upper_bound);
            new_filtered_A1 = filtered_A1(x(i) - Lower_bound:x(i) + Upper_bound);
        end
       
        
        [filtered_r,filtered_lags] = xcorr(new_filtered_A0,new_filtered_A1,MaxLags);
        index_filtered_interp1 = find(filtered_r == max(filtered_r(fix(length(filtered_r)/2:end))),1);      % only looking at correlations that result in (+) wave speeds (1/9/24)
        filtered_difference(i) = filtered_lags(index_filtered_interp1);
        filtered_interp1_speed(i) = ((.01)/(filtered_difference(i)/TimeToIndexConversion));                 % 'abs' was removed (1/9/24)
        if filtered_interp1_speed(i) == inf
            filtered_interp1_speed(i) = p.Results.MinSpeed;
        end
    
        time(i) = x(i)/TimeToIndexConversion;
    end
    
    %disp(filtered_interp1_speed)
    %% removing outliers from wavespeed
    wavespeed_NoOutliers = zeros(1,length(filtered_interp1_speed));              % create a new vector for the wavespeed squared w/o outliers
    time_NoOutliers = zeros(1,length(time));                                        % create a new vector for time to match the wavespeed w/o outliers
    for i = 1:length(filtered_interp1_speed)                                             % walk through the WS^2 vectors
        if filtered_interp1_speed(i) > 200 || filtered_interp1_speed(i) < -200             % all values >50^2 and <10^2 set to NaN
            wavespeed_NoOutliers(i) = NaN;
            time_NoOutliers(i) = NaN;
        else                                                                        % Otherwise, no change
            wavespeed_NoOutliers(i) = filtered_interp1_speed(i);
            time_NoOutliers(i) = time(i);
        end
    end
    wavespeed_NoOutliers = rmmissing(wavespeed_NoOutliers);         % remove the NaNs
    time_NoOutliers = rmmissing(time_NoOutliers);
    
    %max(wavespeed_NoOutliers)
    AvgSpeed = sum(wavespeed_NoOutliers)/length(wavespeed_NoOutliers);
    assignin('base', 'sampleAvgSpeed', AvgSpeed);
    
    
    %% plotting the speed of the signals
    subplot(5,5,n)
    hold on
    plot(time_NoOutliers, wavespeed_NoOutliers,"bx")                                  % plot the wave speeds
    title(filename);
    result = sum(wavespeed_NoOutliers)/length(wavespeed_NoOutliers);
    xlabel('time (s)');
    ylabel('speed (m/s)');
    %ylim([0,60])
    
    hold off
    n = n+1;
end

end
