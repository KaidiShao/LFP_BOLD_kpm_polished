function y = filterSignal_mirror(x, Fs, lo, up, mirror_flag)
% Filter raw signal
%   y = filterSignal(x, Fs) filters the signal x. Each column in x is one
%   recording channel. Fs is the sampling frequency. The filter delay is
%   compensated in the output y.

% design a Butterworth IIR filter with the order 20
hd = design(fdesign.bandpass('N,F3dB1,F3dB2',20,lo,up,Fs),'butter');
% bidirectional filtering to compensate phase delay

if ~mirror_flag
    x2=flipud(filter(hd,x)); % forward filtering
    y=flipud(filter(hd,x2)); % reward filtering
else
    x_origin = x;
    mirror_length = 5*10^5;
    x_mirror = [x(mirror_length:-1:1,:); x; x(end:-1:end-mirror_length+1, :)];

    x2=flipud(filter(hd,x_mirror)); % forward filtering
    y_mirror=flipud(filter(hd,x2)); % reward filtering

    y = y_mirror([1:size(x, 1)]+mirror_length,:);

end


