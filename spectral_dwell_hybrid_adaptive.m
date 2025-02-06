%% Load Data
% Designate paths for audio files and spectral CSV
%audioFilePathStart = "D:\Strobe_Spectra\copper_bell_A.mp3";
%audioFilePathMain = "D:\Strobe_Spectra\DM.mp3";
data = readtable("C:\Users\errat\Desktop\phd\DM\DM_mapping.csv");

% Add required paths for SCCS strobe device interface
addpath('D:\Strobe_Spectra')

%% User Input: Song Name & Dwell Time
song_name = input('What song are you playing? ', 's');

dwell_time = input('Enter dwell time (in seconds): ');
if dwell_time < 1
    dwell_time = 1;
    disp('Dwell time too low, setting to minimum (1 second).');
end

%% User Input: Interval Selection
direction_choice = input('Choose interval direction: (1) Up, (2) Down: ');
intervals = {'Unison', 'Minor Second', 'Major Second', 'Minor Third', 'Major Third', 'Perfect Fourth', 'Tritone', 'Perfect Fifth', 'Minor Sixth', 'Major Sixth', 'Minor Seventh', 'Major Seventh', 'Octave'};

for i = 1:length(intervals)
    fprintf('(%d) %s\n', i, intervals{i});
end
interval_choice = input('Choose interval (1-13): ');

if interval_choice == 1
    column_name = 'Adjusted_Corr_Freq';
else
    if direction_choice == 1
        column_name = sprintf('%s_Up_Freq', strrep(intervals{interval_choice}, ' ', ''));
    else
        column_name = sprintf('%s_Down_Freq', strrep(intervals{interval_choice}, ' ', ''));
    end
end

%% Extract Data & Compute Dwell Time Grouping
time = data.Time;
frequencyValues = data.(column_name);
amplitudeValues = data.Amplitude;
brightness = data.Scaled_Amplitude;

song_length = max(time) - min(time);
frameDurationS = 1/2000;
sampleTimes = (0:frameDurationS:song_length-frameDurationS)';
sampleTimesDwell = floor(sampleTimes / dwell_time) * dwell_time;

uniqueDwellPeriods = unique(sampleTimesDwell);
dominantFreqs = zeros(size(uniqueDwellPeriods));

%% Adaptive Threshold Calculation
window_size = round(length(frequencyValues) * 0.05);
rolling_std = movstd(frequencyValues, window_size);
adaptive_thresholds = rolling_std * 0.5 + movmean(abs(diff([0; frequencyValues])), window_size) * 0.1;

for i = 1:length(uniqueDwellPeriods)
    dwellStart = uniqueDwellPeriods(i);
    dwellEnd = dwellStart + dwell_time;
    inDwellWindow = (time >= dwellStart) & (time < dwellEnd);
    
    if any(inDwellWindow)
        freqs = frequencyValues(inDwellWindow);
        amps = amplitudeValues(inDwellWindow);
        
        weightedFreq = sum(freqs .* amps) / sum(amps);
        modeFreq = mode(freqs);
        
        adaptive_thresh = adaptive_thresholds(i);
        if abs(weightedFreq - modeFreq) < adaptive_thresh
            dominantFreqs(i) = weightedFreq;
        else
            dominantFreqs(i) = max(freqs);
        end
    else
        if i > 1
            dominantFreqs(i) = dominantFreqs(i-1);
        else
            dominantFreqs(i) = 0;
        end
    end
end

interpolatedFreqs = interp1(uniqueDwellPeriods, dominantFreqs, sampleTimesDwell, 'nearest', 'extrap');
interpolatedRawFreqs = interp1(time, frequencyValues, sampleTimes, 'nearest', 'extrap');
brightnessInterpMethod = 'linear';
interpolatedBrightness = round(interp1(time, brightness, sampleTimes, brightnessInterpMethod));

dynamicDutyCycle = (interpolatedBrightness / 255) * 100;
strobe = zeros(size(sampleTimes));
for i = 1:length(sampleTimes)
    strobe(i) = (1 + square(2 * pi * mean(interpolatedFreqs(1:i)) * sampleTimes(i), dynamicDutyCycle(i))) / 2;
end

ledONBitmap = binary8ToUint8(repmat(strobe, 1, 8));
centralBrightness = round(interp1(sampleTimes, interpolatedBrightness, sampleTimes, brightnessInterpMethod));
dacChannelValuesPerSample = zeros(length(sampleTimes), 5);
dacChannelValuesPerSample(:, 1) = centralBrightness;
dacChannelValuesPerSample(:, 2:5) = repmat(centralBrightness, 1, 4);

preparedStrobeData2D = [ledONBitmap, dacChannelValuesPerSample];
preparedStrobeData1D = reshape(preparedStrobeData2D', [], 1)';

comPort = 'COM5';  
filename = [song_name, '.txt'];
success = SCCS_strobe_load_device(preparedStrobeData1D, comPort, filename);

if success
    disp('Strobe data successfully loaded and played on the device.');
else
    disp('Failed to load strobe data to the device.');
end

%% Visualization
figure;
tiledlayout(2,1)

% 🎵 Adaptive Frequency + Brightness
nexttile
title({['Selected Audiovisual Interval (' column_name ')'], ['Dwell Time: ' num2str(dwell_time) 's | Fully Adaptive']})

yyaxis left
plot(sampleTimes, interpolatedFreqs, 'b', 'LineWidth', 1.5, 'DisplayName', 'Adaptive Frequency')
hold on
plot(sampleTimes, interpolatedRawFreqs, 'c--', 'LineWidth', 1, 'DisplayName', 'Raw Frequency')
ylabel("Frequency (Hz)")
xlabel("Time (seconds)")
legend('Location', 'best')

yyaxis right
plot(sampleTimes, dynamicDutyCycle, 'r', 'LineWidth', 1.5, 'DisplayName', 'Brightness / Duty Cycle')
ylabel("Strobe Light Intensity / Duty Cycle (%)")
ylim([0, 100])
legend('Location', 'best')

% ⚡ Strobe Output
nexttile
title('Strobe Output Over Time')
hold on
plot(sampleTimes, strobe, 'k', 'LineWidth', 1.5, 'DisplayName', 'Strobe')
yticks([0, 1])
yticklabels(["Off", "On"])
xlabel("Time (seconds)")
legend('Location', 'best')

hold off;

function value = binary8ToUint8(bitArray)
    value = sum([2^7 2^6, 2^5, 2^4, 2^3, 2^2, 2^1, 2^0] .* bitArray, 2);
end
