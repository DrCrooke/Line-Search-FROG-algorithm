function [Obj, retrieved_trace, errors, BestIteration, BestSigma, Perturbations] = LSF_retrieve(trace,F, D, varargin)
%LSF_retrieve retrieves SHG FROG measurements using the Line-Search
%Retrieval method including sigma check to ensure convergence. Replace the
%error function to change FROG geometry or pulse measurement technique
%(d-scan for example). The algorithm looks for a complex vector with the
%lowest cost (error) function, so by replacing it it will work for other
%geometries.
%
%
% Code author: Christoffer Oxelmark Krook. 
%
%   Returns    :
%       Obj             : Retrieved electric field              - 1xN complex double
%       retrieved_trace : Retrieve trace                        - NxN double
%       errors          : Retrieval errors per iteraiton        - 1xIterations
%       BestIteration   : Iteration with best retrieved field   - numeric
%       BestSigma       : Sigma from best retrieval             - numeric
%       Perturbations   : Number of times sigma check perturbed - numeric
%
%   Arguments  :
%       trace      : Measurement data                   - NxN double         - Required
%       F          : Frequency vector of retrieved field- 1xN double         - Required 
%                                                                              (Shift around zero i.e feed it
%                                                                              fftshift(F) where F goes
%                                                                              from -bandwidth/2 to bandwidth/2)
%       D          : Vector containing delay-values     - 1xN double         - Required
%       Iterations : Number of iterations to perform    - numeric            - optional (2000 default)
%       Verbose    : Verbose optimization               - boolean            - optional (true default)
%       Movie      : Displays retrieval progress        - boolean            - optional (false default)
%       Parfind    : parallel optimization              - boolean            - optional (false default)
%       Lines      : Number of parallel lines to cast   - numeric            - optional (8 default)
%       Guess      : Initial guess                      - 1xN complex double - optional (default is
%                                                                              intensity margin or random field)
%       Sigma      : Sigma check                        - boolean            - optional (true default)

%   Retrieval stops when number of iterations are up or improvement of
%   normalized FROG error over 5 iterations is less than 1e-6.

% assuming NxN-trace, size of retrieved field.
N = length(trace);

% colormap - change to ones preference
colorm = "jet";

% parser for input arguments
p = inputParser;

% define arguments - standard values may be altered here.
addRequired(p, 'trace');
addParameter(p, 'Iterations', 2000);
addParameter(p, 'Verbose', true);
addParameter(p, 'Movie', false);
addParameter(p, 'Parfind', false);
addParameter(p, 'Lines', 8);
addParameter(p, 'Guess', sum(trace, 1))
% addParameter(p, 'Guess', rand(1,N) .* exp(1i*2*pi*rand(1,N))); % initial
% guess is not super important for convergence likelihood,
% but more important for convergence speed.
addParameter(p, 'Sigma', true);

% parse input data
parse(p, trace, varargin{:});
trace          = p.Results.trace;
iterations     = p.Results.Iterations;
verbose        = p.Results.Verbose;
movie          = p.Results.Movie;
parfind        = p.Results.Parfind;
lines          = p.Results.Lines;
x              = p.Results.Guess;
do_sigma_check = p.Results.Sigma;

% check so that start-guess is not outside of permitted parameter space
if max(abs(x)) > 1
    x = x ./max(abs(x)) .* rand; % first scales to unity, and then to a uniformly random number between 0 and 1; 
end

%
%
%
%
% threshold levels. Needs to be tuned to ones preference. 
check_sigma_threshold = 5e-5;
perturbation_threshold = 1e-5;
stop_threshold = 1e-6;
%
%
%
%


% naturally, this is currently the best guess.
bestGuess = x;

% parameter to contain slopes
slopes = zeros(1, iterations/5);

% start iteration
curr_it = 1;

% preallocating errors vector
errors = zeros(1, iterations);

% creates 2D-gaussian filter for sigma check. Parameters may be
% experimented with.
[X, Y] = meshgrid(-N/2:N/2-1,-N/2:N/2-1); % Create coordinates
G = exp(-2.77.*((Y/ceil(N/3)).^2 + (X/ceil(N/3)).^2)); % Create the 2D-Gaussian
G(G < 0.5) = 0; % Everything outside the FWHM is filtered away.

% Are we checking sigma or not? (conditions to start sigma check is
% indication of stagnation. We start with false.)
checkingSigma = false;

% preallocating best error
best_err = error(x, trace, F, D);

% preallocate bestsigma
BestSigma = 100;
sig = NaN;

% preallocate bestIteration
BestIteration = 0;

% preallocating Perturbations
Perturbations = 0;

while curr_it <= iterations
    % Subroutine. This can only be ran after the first iterations
    if curr_it > 1
        x = heuristic(x, xm, xr, xl, xe, xme, xre, xle);
    end

    % perform an iteration
    [xl, xm, x, xr, xle, xme, xe, xre] = parallel_line_search_optimization(x, F, D, lines, trace, N, parfind);
    errors(curr_it) = error(x, trace, F, D); % store retrieved error

    % check if returned vector is best so far, if so, replace the values
    % currently held.
    if errors(curr_it) < best_err
        bestGuess = x;
        best_err = errors(curr_it);
        BestIteration = curr_it;
    end
    % every 5th iteration perform some diagnostic
    if mod(curr_it, 5) == 0
        % check moving error over 5 (slope)
        k = (errors(curr_it) - errors(curr_it - 4)) / (curr_it - (curr_it - 4));
        slopes(curr_it/5) = k;

        % check if retrieval has stagnated
        if abs(slopes(curr_it/5)) < check_sigma_threshold
            checkingSigma = true;
        end

        % perform sigma check
        if checkingSigma & do_sigma_check
            [sig, lowpass_psi] = checkSigma(trace, x, F, D, G); % calculate sigma
            if sig < BestSigma % if sigma is best so far, replace that value
                BestSigma = sig; % Currently the algorithm prioritizes frog error G'
                % over sigma, therefore only this is replaced and not the
                % other "best"-variables.
            end
        end

        % if slope is less than 1e-6, perturb it
        % we must also do some iterations first so that we don't return too
        % soon. This was a problem long ago and I don't remember why, but
        % lets keep it for the sake of nothing. Anyway these checks would
        % not be done before 49 iterations.
        if curr_it > 49
            if abs(slopes(curr_it/5)) < perturbation_threshold && BestSigma > 5 && do_sigma_check
                if mod(Perturbations, 3) == 0
                    x = rand.*exp(-2.77.*((-N/2:N/2-1)./(abs(N*randn*0.3))).^2).*exp(1i.*unwrap(angle(x))); % random gaussian intensity and keep the current phase
                else
                    x = rand.*exp(-2.77.*((-N/2:N/2-1)./(abs(N*randn*0.3))).^2).*exp(1i.*randn(1, N)); % start over from scratch.
                end
                Perturbations = Perturbations + 1;
                checkingSigma = false;
                if max(abs(x)) > 1
                    disp("")
                end
            elseif abs(slopes(curr_it/5)) < stop_threshold % if other criterias are met, we are
                % done and therefore we return early before using all iterations.
                [sig, lowpass_psi] = checkSigma(trace, x, F, D, G);
                break;
            end
        end
    end

    %
    if mod(curr_it, 25) == 0 && curr_it > 5
        if verbose % write some stuff every 25th iteration
            clc
            disp("Error         : " + num2str(errors(curr_it)));
            disp("Iteration     : " + num2str(curr_it));
            disp("Slope         : " + num2str(slopes(curr_it/5)));
            disp("checkingSigma : " + checkingSigma + ", Sigma = " + num2str(sig))
            disp("Best it       : " + num2str(BestIteration) + ", Best error: " + num2str(best_err) + ", Best Sigma: " + num2str(BestSigma));
        end
    end
    if mod(curr_it, 25) == 0 && curr_it > 5
        if movie % display some stuff every 25th iteration
            temptrace = trace;
            figure(5); clf; set(gcf, 'Position',  [100, 100, 1000, 700])
            if checkingSigma && do_sigma_check
                subplot(2,4,1); imagesc(D.*1e15, fftshift(F).*1e-12, temptrace); title("Measurement");colormap(colorm); xlabel("Delay (fs)"); ylabel("Normalized frequency (THz)");
                subplot(2,4,2); imagesc(D.*1e15, fftshift(F).*1e-12, createProperTrace(x, F, D)); title("Retrieved trace");colormap(colorm); xlabel("Delay (fs)"); ylabel("Normalized frequency (THz)");
                subplot(2,4,3); imagesc(D.*1e15, fftshift(F).*1e-12, createProperTrace(bestGuess, F, D)); title("Best retrieval trace");colormap(colorm); xlabel("Delay (fs)"); ylabel("Normalized frequency (THz)");
                subplot(2,4,4); imagesc(D.*1e15, fftshift(F).*1e-12, abs(lowpass_psi)); colorbar; colormap(colorm); title("Current sigma " + num2str(sig)); xlabel("Delay (fs)"); ylabel("Normalized frequency (THz)");
            else
                subplot(2,3,1); imagesc(D.*1e15, fftshift(F).*1e-12, temptrace); title("Measurement");colormap(colorm); xlabel("Delay (fs)"); ylabel("Normalized frequency (THz)");
                subplot(2,3,2); imagesc(D.*1e15, fftshift(F).*1e-12, createProperTrace(x, F, D)); title("Retrieved trace");colormap(colorm); xlabel("Delay (fs)"); ylabel("Normalized frequency (THz)");
                subplot(2,3,3); imagesc(D.*1e15, fftshift(F).*1e-12, createProperTrace(bestGuess, F, D)); title("Best retrieval trace");colormap(colorm); xlabel("Delay (fs)"); ylabel("Normalized frequency (THz)");
            end
            subplot(4,1,3); semilogy(errors); hold on; ylabel("G'"); xlabel("Iteration #"); title("Cost function");% hold on; plot((0:1:iterations), k*(0:1:iterations) + intercept); ylim([0 1]); xlim([0 curr_it]);
            subplot(4,1,4); yyaxis left; plot(D.*1e15, abs(bestGuess).^2, 'Color', 'red', 'LineWidth', 1.5); xlabel("Time"); ylabel("Intensity (A.U)");
            hold on; plot(D.*1e15, abs(x).^2, 'Color', 'black', 'LineWidth', 1.5); hold off;
            yyaxis right; plot(D.*1e15, unwrap(angle(x)), '--'); title("Retrieved pulse intensity and phase"); ylabel("Temporal phase (rad)"); yyaxis left; grid on; set(gca, 'LineWidth', 2)
            legend('Best retrieval', 'Current retrieval', 'Current retrieval phase');
        end
    end

    

    curr_it = curr_it + 1;
end
if do_sigma_check == false
    BestSigma = sig;
end
Obj = bestGuess./max(abs(bestGuess));
retrieved_trace = createProperTrace(Obj, F, D);
errors = errors(errors ~= 0);
clc;
disp("Retrieval finished")
end

% cost function
function err = error(x, trace, F, D)
retrieved_trace = createProperTrace(x, F, D); % generate a trace based on current guess
mu = sum(sum(trace.*retrieved_trace))./sum(sum(retrieved_trace.^2)); % re-scaling factor
err = sqrt(sum(sum((trace - mu.*retrieved_trace).^2)))./sqrt(sum(sum(trace.^2))); % calculate the normalized frog error
end

% subroutine to check if the function could not find a better vector
function x = heuristic(x, xm, xr, xl, xe, xme, xre, xle)
if min([xe xme xre xle]) ~= xe % check if returned error was not the lowest.
    % if not, return the one with the 2nd lowest error of the remaining
    % vectors.
    if min([xme xre xle]) == xle
        x = xl;
    elseif min([xme xre xle]) == xme
        x = xm;
    else
        x = xr;
    end
end
end

% generate lines
function[xl, xm, xr] = selectPoints(x, N)
% this method is currently a mess, BUT gets the job done. Refinement
% possible. Also bruteforce by incrementally increasing ss_1 and ss_2 until threshold is an
% alternative, however without elegance.

d1 = (sign(randn(1,N)).*rand(1,N) + 1i.*sign(randn(1,N)).*rand(1,N)); % random directional vector
d1 = d1./(abs(d1));% normalize.
% d = rand(1,N) .* exp(1i*2*pi*rand(1,N));

% find values of s so that one element in xl and xr, strictly, has a
% magnitude of unity
a = real(x);
b = imag(x);
c = real(d1);
d = imag(d1);

s1 = min((-(2.*a.*c + 2.*b.*d) + sqrt((2.*a.*c+2.*b.*d).^2 - 4.*(c.^2+d.^2).*(a.^2+b.^2-1)))./(2.*(c.^2+d.^2)));
s2 = -min(abs(((-(2.*a.*c + 2.*b.*d) - sqrt((2.*a.*c+2.*b.*d).^2 - 4.*(c.^2+d.^2).*(a.^2+b.^2-1)))./(2.*(c.^2+d.^2)))));

xl = x + s1.*d1;
xr = x + s2.*d1;
xm = (xl+xr)./2;
end

function [sig, fftransformed] =  checkSigma(trace, x, F, D, G)
% low pass filter, and compare noise to peak structure present.
N = length(x);
current_guess_trace = createProperTrace(x, F, D);
mu = sum(sum(trace.*current_guess_trace))./sum(sum(current_guess_trace.^2));
tracediff = trace - mu.*current_guess_trace; % psi
ffttrace = fftshift((fft2(tracediff)));
fftransformed = (ifft2(fftshift(G.*ffttrace))); % psi_LP.

floorval = (std2(fftransformed(1:floor(N/4), 1:floor(N/4))) + std2(fftransformed(1:floor(N/4), end-(floor(N/4)-1):end)) + std2(fftransformed(end-(floor(N/4)-1):end, 1:floor(N/4))) + std2(fftransformed(end-(floor(N/4)-1):end, end-(floor(N/4)-1):end)))/4;
sig =  max(max(abs(fftransformed)))/floorval;
end



% Line-search loop
function [x, xl, xm, xr, xe, xme, xle, xre] = falsifyLine(xl, xm, xr, x, trace, F, D)
iterations_without_improvement = 0;

% initial error values for the fields we start with
xle = error(xl, trace, F, D);
xre = error(xr, trace, F, D);
xme = error(xm, trace, F, D);
xe = error(x, trace, F, D);


% how many times should be reduce the line? This gives a random value
% between 1 and 12 in an attempt to combine a rough check (few reductions)
% and a detailed check (many reductions) in an attempt to cover as much of
% the parameter space as possible. Should reduce local minima sensitivity
% slightly.
reductions = ceil(12*rand);


% main function of line-search. Should be as https://research.chalmers.se/publication/526076/file/526076_Fulltext.pdf
% but with added on-the-spot error calculations to reduce total number of
% calls to error().
while iterations_without_improvement < reductions
    if xle < xre && xle < xme
        x_new = (xl + xm)./2;
        xne = error(x_new, trace, F, D);
        if xne < xle
            x = x_new;
            xe = xne;
        end
        xr = xm;
        xre = xme;

        xm = x_new;
        xme = xne;
    elseif xre < xle && xre < xme
        x_new = (xm + xr)./2;
        xne = error(x_new, trace, F, D);
        if xne < xre
            x = x_new;
            xe = xne;
        end
        xl = xm;
        xle = xme;
        xm = x_new;
        xme = xne;
    else
        if sqrt(sum((xl-xm).^2)) >= sqrt(sum((xm-xr).^2))
            x_new = (xl + xm)./2;
            xne = error(x_new, trace, F, D);
            if xne < xme
                x = x_new;
                xe = xne;
                xr = xm;
                xre = xme;
                xm = x_new;
                xme = xne;
            else
                xl = x_new;
                xle = xne;
            end
        else
            x_new = (xm+xr)./2;
            xne = error(x_new, trace, F, D);
            if xne < xme
                x = x_new;
                xe = xne;
                xl = xm;
                xle = xme;
                xm = x_new;
                xme = xne;
            else
                xr = x_new;
                xre = xne;
            end
        end
    end

    % If we reduced down to the same point.
    if xle == xme || xre == xme
        return
    end

    % better point?
    if x == x_new
        iterations_without_improvement = 0;
    else
        iterations_without_improvement = iterations_without_improvement + 1;
    end
end

errs = [xe, xne, xle, xme, xre]; %values are given as output

% sets the best vector found as x
if min(errs) == xle
    x = xl;
    xe = xle;
elseif min(errs) == xme
    x = xm;
    xe = xme;
elseif min(errs) == xr
    x = xr;
    xe = xre;
elseif min(errs) == xne % xne = x_new_error.
    x = x_new;
    xe = xne;
end
end


% this function performs the proper function calls and handles if we want
% to do it in parallel or not.
function [xl, xm, x, xr, xle, xme, xe, xre] = parallel_line_search_optimization(x, F, D, lines, trace, N, parfind)
% Initialize storage for the results
xls = zeros(N, lines);
xms = zeros(N, lines);
xrs = zeros(N, lines);
xes = zeros(1, lines);
xles = zeros(1, lines);
xres = zeros(1, lines);
xmes = zeros(1, lines);
xtemps = zeros(N, lines);
% Perform the line-search optimization in parallel
if parfind
    parfor i = 1:lines
        % Generate points
        [xl, xm, xr] = selectPoints(x, N);

        % Perform direct search
        [tempx, xl1, xm1, xr1, xe, xme, xle, xre] = falsifyLine(xl, xm, xr, x, trace, F, D);

        % Store results
        xls(:, i) = xl1;
        xms(:, i) = xm1;
        xrs(:, i) = xr1;
        xtemps(:, i) = tempx;
        xes(i) = xe;
        xles(i) = xle;
        xres(i) = xre;
        xmes(i) = xme;
    end
else
    % same as above, however with non-parallell for loop.
    for i = 1:lines
        % Generate points
        [xl, xm, xr] = selectPoints(x, N);
        % plotLines([xl(1) xl(50) xl(78)], [xm(1) xm(50) xm(78)], [x(1) x(50) x(78)], [xr(1) xr(50) xr(78)]);
        % Perform direct search
        [tempx, xl1, xm1, xr1, xe, xme, xle, xre] = falsifyLine(xl, xm, xr, x, trace, F, D);
        % Store results
        xls(:, i) = xl1;
        xms(:, i) = xm1;
        xrs(:, i) = xr1;
        xtemps(:, i) = tempx;
        xes(i) = xe;
        xles(i) = xle;
        xres(i) = xre;
        xmes(i) = xme;
    end
end

% Find the best result with the lowest xe
[~, idx] = min(xes);

% Select the best results
xl = xls(:, idx).';
xm = xms(:, idx).';
xr = xrs(:, idx).';
x = xtemps(:, idx).';
xe = xes(idx);
xle = xles(idx);
xre = xres(idx);
xme = xmes(idx);
end
