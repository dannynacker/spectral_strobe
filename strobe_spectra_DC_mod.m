% ------------------------------------------------------------------------------ 
%  This MATLAB script is designed to be used with the spectral.ipynb process 
%  which pulls spectral data from an uploaded audio file and converts dominant 
%  frequency as determined by amplitude per time slice into their corresponding 
%  "notes" in a different octave to be used in a stroboscopic experience.
%  Stroboscopic frequency and amplitude per time slice remains a function of those 
%  of the auditory data. These data are then fed to the SCCS Strobe Device.
%  Upon running this script, you will be prompted to select a shift in the
%  interval relationship between the photic and auditory stimuli. 
%  Before the file is sent to the light, brightness, corresponding frequency, 
%  our selected transposed relative frequency, and 
%  their associated intermodulation distortion components are plotted 
%  alongside the square waveform output of the light. Once the file has been read 
%  by the light, a bell tone will play to signal to the participant to shut their 
%  eyes, relax, and enjoy the show!
%  This script differs from strobe_spectra_intervals because now the
%  DUTY CYCLE of the stroboscope varies as a function of amplitude, 
%  creating a further modulation of subjective intensity that anecdotal reports
%  show a greater preference for compared to that of a constant 50% duty cycle
% ------------------------------------------------------------------------------ 

% Designate the paths to your onset signal and music files 
audioFilePathStart = "D:\Strobe_Spectra\copper_bell_A.mp3";
audioFilePathMain = "D:\Strobe_Spectra\jupiter.mp3";

% Designate the path to your spectral CSV file
data = readtable("D:\Strobe_Spectra\luca\luther.csv");
             'Minor Seventh', 'Major Seventh', 'Octave'};

% Adjust the loop to start from 1 instead of 0
for i = 1:length(intervals)
    fprintf('(%d) %s\n', i, intervals{i});
end
interval_choice = input('Choose interval (1-13): ');

% Construct the column name based on the user input
if interval_choice == 1
    column_name = 'Adjusted_Corr_Freq';
    disp('Using Unison (Corresponding Frequencies).');
else
    if direction_choice == 1
        column_name = sprintf('%s_Up_Freq', strrep(intervals{interval_choice}, ' ', ''));
        disp(['Using ', intervals{interval_choice}, ' Up.']);
    elseif direction_choice == 2
        column_name = sprintf('%s_Down_Freq', strrep(intervals{interval_choice}, ' ', ''));
        disp(['Using ', intervals{interval_choice}, ' Down.']);
    else
        error('Invalid direction choice! Please choose 1 or 2.');
    end
end

% Ensure the constructed column name is correct
fprintf('Accessing column: %s\n', column_name);

% Extract the data based on the selected column
time = data.Time;
frequencyValues = data.(column_name);
adjustedCorrFreq = data.Adjusted_Corr_Freq;  % The corresponding frequency to compare against
brightness = data.Scaled_Amplitude;

% Calculate the song length in seconds
song_length = max(time) - min(time);

% Print the song length
disp(['Song length: ', num2str(song_length), ' seconds']);

% Generate sample times
frameDurationS = (1/2000); % Time duration of each frame
sampleTimes = (0:frameDurationS:song_length-frameDurationS)'; % Generate a list of sample timestamps

% Interpolate the data based on sample times
frequencyInterpMethod = 'nearest'; 
brightnessInterpMethod = 'linear';

interpolatedFreqs = interp1(time, frequencyValues, sampleTimes, frequencyInterpMethod);
interpolatedCorrFreqs = interp1(time, adjustedCorrFreq, sampleTimes, frequencyInterpMethod);
interpolatedBrightness = round(interp1(time, brightness, sampleTimes, brightnessInterpMethod));

% Generate strobe waveform
avgFreqSinceStart = cumsum(interpolatedFreqs) ./ (1:length(interpolatedFreqs))';

%% Step 4: Generate strobe waveform
% Convert scaled amplitude (0-255) to duty cycle (0-100%)
dynamicDutyCycle = (interpolatedBrightness / 255) * 100;

% Preallocate the strobe waveform
strobe = zeros(size(sampleTimes));

% Generate the strobe waveform dynamically
for i = 1:length(sampleTimes)
    strobe(i) = (1 + square(2 * pi * avgFreqSinceStart(i) * sampleTimes(i), dynamicDutyCycle(i))) / 2;
end

%% Step 5: Create LED ON bitmap
ledONBitmap = binary8ToUint8(repmat(strobe, 1, 8)); % 8-bit bitmap

% Step 6: Interpolate central and ring brightness
centralBrightness = round(interp1(sampleTimes, interpolatedBrightness, sampleTimes, brightnessInterpMethod));
ringBrightness = round(interp1(sampleTimes, interpolatedBrightness, sampleTimes, brightnessInterpMethod));

% Step 7: Preallocate DAC channel values
dacChannelValuesPerSample = zeros(length(sampleTimes), 5);
dacChannelValuesPerSample(:, 1) = centralBrightness;
dacChannelValuesPerSample(:, 2:5) = repmat(ringBrightness, 1, 4);

% Step 8: Generate 2D and 1D arrays
preparedStrobeData2D = [ledONBitmap, dacChannelValuesPerSample];
preparedStrobeData1D = reshape(preparedStrobeData2D', [size(preparedStrobeData2D, 1) * size(preparedStrobeData2D, 2), 1])';

% preparedStrobeData1D now contains the final 1D array.

% Calculate intermodulation distortion components
differenceTones = abs(interpolatedFreqs - interpolatedCorrFreqs);
summationTones = interpolatedFreqs + interpolatedCorrFreqs;

% Clean up the interval and direction names for the title and legend
intervalNameParts = strsplit(column_name, '_');
intervalName = intervalNameParts{1}; % Extract the interval name

if length(intervalNameParts) > 2
    direction = intervalNameParts{2}; % Extract the direction (Up or Down)
    directionFull = [upper(direction(1)), direction(2:end)]; % Capitalize the direction
    clean_column_name = [intervalName, ' ', directionFull]; % Combine interval and direction
else
    clean_column_name = intervalName; % For the unison case, just use the interval name
end

% Adjust for the case where the interval is "Unison"
if strcmp(intervalName, 'Unison')
    clean_column_name = 'Unison';
end

% Compute the maximum y value from all left-axis plotted data
maxY = max([max(interpolatedFreqs), max(interpolatedCorrFreqs), max(differenceTones), max(summationTones)]);

% Plot the results
figure;
tiledlayout(2,1)

%% Plot the selected frequency and IMD Components
nexttile
title(['Selected Audiovisual Interval (' clean_column_name ') and IMD Components'])
yyaxis left
plot(sampleTimes, interpolatedFreqs, 'b', 'DisplayName', clean_column_name)
hold on
plot(sampleTimes, interpolatedCorrFreqs, 'c--', 'DisplayName', 'Music Frequency')
plot(sampleTimes, differenceTones, 'g--', 'DisplayName', 'Difference Tone')
plot(sampleTimes, summationTones, 'm--', 'DisplayName', 'Summation Tone')
ylabel("Frequency (Hz)")
ylim([0, maxY])
xlabel("Time (seconds)")
xlim([0, song_length])
legend('Location', 'best')

yyaxis right
plot(sampleTimes, dynamicDutyCycle, 'r', 'DisplayName', 'Duty Cycle & Brightness (%)') % <-- Added Duty Cycle Plot
ylabel("Strobe Light Intensity / Duty Cycle (%)")
ylim([0, 100])
legend('Location', 'best')

%% Plot Strobe Output
nexttile
title('Strobe Output')
hold on
plot(sampleTimes, strobe, 'DisplayName', 'Strobe')
ylim([-0.5,1.5])
yticks([0,1]);
yticklabels(["Off", "On"])
xlabel("Time (seconds)")
xlim([0, song_length])
legend('Location', 'best')

pause(1);

% Specify COM port and filename
comPort = 'COM5'; % Replace with your actual COM port
filename = [song_name, '.txt']; % Ensure this matches the device's 8.3 filename requirement

% Load the strobe data to the device
success = SCCS_strobe_load_device(preparedStrobeData1D, comPort, filename);

if success
    % Play the main audio file
    disp('Strobe data successfully loaded and played on the device.');
else
    disp('Failed to load strobe data to the device.');
end

function value = binary8ToUint8(bitArray)
    value = sum([2^7 2^6, 2^5, 2^4, 2^3, 2^2, 2^1, 2^0] .* bitArray, 2);
    return;
end
