% Code author: Christoffer Oxelmark Krook. 19/5-2025
close all
% Top-level file to run the LSF algorithm.
% The Sigma-Check invokes the std2- function of the Image Processing
% Toolbox. To use proper parallell execution the paralell processing
% toolbox also needs to be installed. The rest should be included already.

% For help how to run LSF the way you'd like, write "help LSF_retrieve" in
% the console
clc; clear all;

% set up grid.
N = 128;
time = linspace(-100e-15, 100e-15, N);
dt = time(2) - time(1);
F = (-N/2:N/2-1).*(1/(N*dt));

% load pulses
pulse_bank = load("pulses.mat");
pulse_number = 1; % chooise value between 1 and 250.
pulse = pulse_bank.pulses(:, pulse_number);

% create trace
trace = createProperTrace(pulse.', fftshift(F), time);
% add noise
noise_amount = 15; % in percent.
trace = trace + randn(N,N).*(noise_amount/100.*max(max(abs(trace))));% LSF paper style noise application


% Retrieval parameters
perform_sigma_check = true; % do you want to run with the sigma check? Convergence improving step
iterations = 50000; % Maximum number of allowed iteratioins
parcheck = true; % check with paralell for loop
lines = 14; % number of paralell lines to evaluate.
movie = true; % display retrieval live
verbose = true; % print to console.

%
%
% Thresholds in line 82-84 in LSF_retrieve might need tuning to ones
% specific case
%
%

% retrieve the pulse with LSF. Write "help LSF_retrieve" in console for
% more detailed information.
[Obj, retrieved_trace, errors, BestIteration, BestSigma, Perturbations] = LSF_retrieve(trace, fftshift(F), time, 'Iterations', iterations, 'Movie', movie, 'Parfind', parcheck, 'Lines', lines, 'Sigma', perform_sigma_check);
clc;
close all;


% remove ambiguities, interpolate fields first so the accuracy is higher
% with ambiguous time delay. Function borrowed from https://doi.org/10.1364/OPTICA.3.001320
pulse = interp1(time, pulse, linspace(-100e-15, 100e-15, 2048), 'pchip', 0);
Obj = interp1(time, Obj, linspace(-100e-15, 100e-15, 2048), 'pchip', 0);

time = linspace(-100e-15, 100e-15, 2048);

disp("Removing ambiguities, please wait ...")
[Obj, pulse] = best_sol(Obj.', pulse.'); % this function corrects time-direction and amboguous time-shift.

% temporal intensities
ref_int = abs(pulse).^2;
ret_int = abs(Obj).^2;

% temporal phase
ref_phase = unwrap(angle(pulse));
ret_phase = unwrap(angle(Obj));

% temporal int error
int_err = sqrt(1/2048.*sum(abs(ref_int - ret_int).^2));
int_err_bins = sqrt(1/(length(ref_int(ref_int > 0.005)))*sum(abs(ref_int(ref_int > 0.005) - ret_int(ref_int > 0.005)).^2));

% spectral intensity
ref_int_spec = abs(fftshift(fft(pulse))).^2; old_val = max(ref_int_spec); ref_int_spec = ref_int_spec./old_val; quote = 1 / old_val;

% rescale them to unity.
ret_int_spec = abs(fftshift(fft(Obj))).^2 .* quote;

int_err_spec = sqrt(1/2048*sum(abs(ref_int_spec - ret_int_spec).^2));
int_err_spec_bins = sqrt(1/length(ref_int_spec(ref_int_spec > 0.005))*sum(abs(ref_int_spec(ref_int_spec > 0.005) - ret_int_spec(ref_int_spec > 0.005)).^2));

% new frequency-vector for plotting purposes, now that the time-vector has
% a higher sampling rate.
dt = time(2) - time(1);
F = (-2048/2:2048/2-1).*(1/(2048*dt));


% plot intensities and phases
figure; subplot(2,3,1); plot(time.*1e15, abs(pulse).^2, 'LineWidth', 2, 'Color', 'black'); title("Temporal intensity"); xlabel("Time (fs)"); ylabel("Normalized intensity (A.U)");
hold on; plot(time.*1e15, abs(Obj).^2, 'LineWidth', 1.5, 'Color', 'red'); ylim([-0.1 1.1]); xlabel("Time (fs)"); ylabel("Normalized intensity (A.U)");
subplot(2,3,2); plot(time.*1e15, ref_phase - ref_phase(1024), 'LineWidth', 2, 'Color', 'black'); ylabel("Phase (rad)");
hold on; plot(time.*1e15, ret_phase - ret_phase(1024), '-', 'LineWidth', 1.5, 'Color', 'red'); ylabel("Phase (rad)"); xlabel("Time (fs)"); title("Temporal phase");
subplot(2,3,3); plot(F.*1e-15, abs(fftshift(fft(pulse))).^2, 'LineWidth', 2, 'Color', 'black'); hold on; title("Spectral intensity"); xlabel("Frequency (PHz)"); ylabel("Intensity (A.U)");
plot(F.*1e-15, abs(fftshift(fft(Obj))).^2, 'LineWidth', 2, 'Color', 'red'); xlim([-0.3 0.3]); ylim([-0 max(abs(fftshift(fft(pulse))).^2).*1.1]); xlabel("Frequency (PHz)"); ylabel("Intensity (A.U)");

subplot(2,3,4); imagesc(time.*1e15, F.*1e-15, trace); title("Input trace"); xlabel("Time (fs)"); ylabel("Frequency (PHz)");
subplot(2,3,5); imagesc(time.*1e15, F.*1e-15, retrieved_trace); title("Retrieved trace"); xlabel("Time (fs)"); ylabel("Frequency (PHz)");
subplot(2,3,6); imagesc(time.*1e15, F.*1e-15, trace - sum(sum(trace.*retrieved_trace))/sum(sum(retrieved_trace.^2)).*retrieved_trace); title("Difference"); xlabel("Time (fs)"); ylabel("Frequency (PHz)");

% angle
delta = acos(abs(dot(Obj, pulse))/sqrt(dot(Obj, Obj)*dot(pulse, pulse)));

% print metrics.
disp(" ")
if BestSigma < 5 && perform_sigma_check == true
    disp("Sigma check indicates convergence to correct solution.")
elseif BestSigma >= 5 && perform_sigma_check == true
    disp("Sigma check indicates convergence to incorrect solution.")
elseif perform_sigma_check == false
    disp("Convergence not evaluated.")
end
disp("Retrieved in " + num2str(length(errors(errors~= 0))) + " iterations")
disp("G' " + num2str(min(errors(errors~=0))))
disp("δ " + num2str(delta) + " radians.");
disp("Temporal ϵ " + num2str(round(int_err.*100, 2)) + " %. Temporal ϵ of relevant bins " + num2str(int_err_bins * 100) + " %")
disp("Spectral ϵ " + num2str(round(int_err_spec.*100, 2)) + " %. Spectral ϵ of relevant bins " + num2str(int_err_spec_bins * 100) + " %");
