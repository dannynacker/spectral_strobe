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
% ------------------------------------------------------------------------------ 

% Designate the paths to your onset signal and music files 
audioFilePathStart = "D:\Strobe_Spectra\copper_bell_A.mp3";
audioFilePathMain = "D:\Strobe_Spectra\jupiter.mp3";

% Designate the path to your spectral CSV file
data = readtable("D:\Strobe_Spectra\jupiter.csv");

% Designate the path to SCCS_strobe_load_device and StrobeDevice.m here for the success function
% and for interfacing with the strobe device
addpath('D:\Strobe_Spectra')

% Ask for the song name
song_name = input('What song are you playing? ', 's');

% Prompt the user to select interval direction
direction_choice = input('Choose interval direction: (1) Up, (2) Down: ');

% Prompt the user to select the interval (semitone shift)
intervals = {'Unison', 'Minor Second', 'Major Second', 'Minor Third', ...
             'Major Third', 'Perfect Fourth', 'Tritone', ...
             'Perfect Fifth', 'Minor Sixth', 'Major Sixth', ...
             'Minor Seventh', 'Major Seventh', 'Octave'};

for i = 0:length(intervals)-1
    fprintf('(%d) %s\n', i, intervals{i+1});
end
interval_choice = input('Choose interval (0-12): ');

% Construct the column name based on the user input
if interval_choice == 0
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
strobeDutyCycle = 50;
strobe = (1 + square(sampleTimes * 2 * pi .* avgFreqSinceStart, strobeDutyCycle)) ./ 2;

% Calculate intermodulation distortion components
differenceTones = abs(interpolatedFreqs - interpolatedCorrFreqs);
summationTones = interpolatedFreqs + interpolatedCorrFreqs;

% Clean up the interval and direction names for the title and legend
clean_column_name = strrep(column_name, '_', ' ');

% Plot the results
figure;
tiledlayout(2,1)

% Plot the selected frequency and IMD Components
nexttile
title(['Selected Audiovisual Interval (' clean_column_name ') and IMD Components'])
yyaxis left
plot(sampleTimes, interpolatedFreqs, 'b', 'DisplayName', clean_column_name)
hold on
plot(sampleTimes, differenceTones, 'g--', 'DisplayName', 'Difference Tone')
plot(sampleTimes, summationTones, 'm--', 'DisplayName', 'Summation Tone')
ylabel("Frequency (Hz)")
xlabel("Time (seconds)")
xlim([0, song_length])
legend('Location', 'best')
hold on
yyaxis right
plot(sampleTimes, interpolatedBrightness, 'color', [1, 0.5, 0], 'DisplayName', 'Brightness')
ylabel("Strobe Light Intensity")
ylim([0,255])

% Plot Strobe Output
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

% Load and play the audio signal to signal the start of the sequence
[yStart, FsStart] = audioread(audioFilePathStart);

pause(1);

% Play the start signal
sound(yStart, FsStart);
pause(5);

% Specify COM port and filename
comPort = 'COM6'; % Replace with your actual COM port
filename = [song_name, '.txt']; % Ensure this matches the device's 8.3 filename requirement

% Load the strobe data to the device
success = SCCS_strobe_load_device(preparedStrobeData1D, comPort, filename);

if success
    % Play the main audio file
    [yMain, FsMain] = audioread(audioFilePathMain);
    sound(yMain, FsMain);

    disp('Strobe data successfully loaded and played on the device.');
else
    disp('Failed to load strobe data to the device.');
end

function value = binary8ToUint8(bitArray)
    value = sum([2^7 2^6, 2^5, 2^4, 2^3, 2^2, 2^1, 2^0] .* bitArray, 2);
    return;
end
