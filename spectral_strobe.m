% ------------------------------------------------------------------------------ 
%  This MATLAB script is designed to be used with the spectral.ipynb process 
%  which pulls spectral data from an uploaded audio file and converts dominant 
%  frequency as determined by amplitude per time slice into their corresponding 
%  "notes" in a different octave to be used in a stroboscopic experience.
%  Stroboscopic frequency and amplitude per time slice remains a function of those 
%  of the auditory data. These data are then fed to the SCCS Strobe Device.
%  Upon running this script, you will be prompted to select either photic 
%  information that corresponds in unison with the auditory information or 
%  frequencies that have been "transposed down a major third" to potentially 
%  generate greater intermodulation distortion components.
%  Before the file is sent to the light, brightness, corresponding frequency, 
%  relative frequency, and intermodulation distortion components are plotted 
%  alongside the square waveform output of the light. Once the file has been read 
%  by the light, a bell tone will play to signal to the participant to shut their 
%  eyes, relax, and enjoy the show!
% ------------------------------------------------------------------------------ 

% Designate the paths to your onset signal and music files 
audioFilePathStart = "D:\new_interface\Experimental_Script\copper_bell_A.mp3";
audioFilePathMain = "D:\spotify_API\jupiter\jupiter.mp3";

% Designate the path to your spectral CSV file
data = readtable("D:\Strobe_Spectra\notated_freq_amp_time_bach.csv");

% Designate the path to SCCS_strobe_load_device here for the success function
addpath('D:\Strobe_Spectra')

% Ask for the song name
song_name = input('What song are you playing? ', 's');

% Prompt the user to choose between corresponding frequencies or relative frequencies
choice = input('Choose frequency type: (1) Corresponding Frequencies, (2) Relative Frequencies: ');

% Extract the data based on the user's choice
time = data.Time;

if choice == 1
    frequencyValues = data.Adjusted_Corr_Freq;
    otherFreq = data.Adjusted_Rel_Freq;
    disp('Using Corresponding Frequencies.');
    freqLabel = 'Corr Freq';
elseif choice == 2
    frequencyValues = data.Adjusted_Rel_Freq;
    otherFreq = data.Adjusted_Corr_Freq;
    disp('Using Relative Frequencies.');
    freqLabel = 'Rel Freq';
else
    error('Invalid choice! Please choose either 1 or 2.');
end

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
interpolatedOtherFreqs = interp1(time, otherFreq, sampleTimes, frequencyInterpMethod);
interpolatedBrightness = round(interp1(time, brightness, sampleTimes, brightnessInterpMethod));

% Generate strobe waveform
avgFreqSinceStart = cumsum(interpolatedFreqs) ./ (1:length(interpolatedFreqs))';
strobeDutyCycle = 50;
strobe = (1 + square(sampleTimes * 2 * pi .* avgFreqSinceStart, strobeDutyCycle)) ./ 2;

% Calculate intermodulation distortion components
differenceTones = abs(interpolatedFreqs - interpolatedOtherFreqs);
summationTones = interpolatedFreqs + interpolatedOtherFreqs;

% Create LED ON/OFF bitmap
ledONBitmap = binary8ToUint8(repmat(strobe, 1, 8));

% Preallocate dacChannelValuesPerSample for efficiency
dacChannelValuesPerSample = zeros(length(sampleTimes), 5);

% Assign brightness values to DAC channels
dacChannelValuesPerSample(:, 1) = interpolatedBrightness;
dacChannelValuesPerSample(:, 2:5) = repmat(interpolatedBrightness, 1, 4);

% Prepare strobe data
preparedStrobeData2D = [ledONBitmap, dacChannelValuesPerSample];
preparedStrobeData1D = reshape(preparedStrobeData2D', [size(preparedStrobeData2D, 1) * size(preparedStrobeData2D, 2), 1])';

% Save the strobe data to a .mat file
save([song_name, '.mat'], 'preparedStrobeData1D');

% Load and play the audio signal to signal the start of the sequence
[yStart, FsStart] = audioread(audioFilePathStart);

disp("Re-reading file list.")
fileListAfterWrite = device.getFileList(); % Refresh the file list after writing

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

% Plot the results
figure;
tiledlayout(2,1)

% Plot the selected frequency and IMD Components
nexttile
title(['Selected Frequency (' freqLabel ') and IMD Components'])
yyaxis left
plot(sampleTimes, interpolatedFreqs, 'b', 'DisplayName', freqLabel)
hold on
plot(sampleTimes, differenceTones, 'g--', 'DisplayName', 'Difference Tone')
plot(sampleTimes, summationTones, 'm--', 'DisplayName', 'Summation Tone')
ylabel("Frequency (Hz)")
xlabel("Time (seconds)")
xlim([0, song_length])
legend show
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
legend show

function value = binary8ToUint8(bitArray)
    value = sum([2^7 2^6, 2^5, 2^4, 2^3, 2^2, 2^1, 2^0] .* bitArray, 2);
    return;
end
