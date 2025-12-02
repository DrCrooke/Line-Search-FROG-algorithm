function trace = createProperTrace3(x, F, D)
%   Standard fourier-domain delay of temporal signal.
N = length(x);
fieldFFT_col = fft(x(:)); % Column vector
D_row = D(:)';                 % Row vector
F_col = F(:);            % Column vector

phase_shifts_matrix = exp(1i * 2 * pi * (F_col * D_row)); 

fieldFFT_matrix = repmat(fieldFFT_col, 1, length(D_row));

P_matrix = ifft(fieldFFT_matrix .* phase_shifts_matrix);

field_col = x(:);
P_times_field = P_matrix .* field_col; % field_col expands across columns

trace = fftshift(abs(fft(P_times_field)/N).^2, 1);
end