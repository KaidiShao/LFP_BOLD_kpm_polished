function [tmpall, freqs, timesout,tmpangle] = timefreqMB(data, srate, varargin)
% TIMEFREQMB - compute time/frequency decomposition of data trials.
%
% Usage:
%     >> [tf, freqs, times]          = timefreq(data, srate);
%     >> [tf, freqs, times, itcvals] = timefreq(data, srate, ...
%                                        'key1', 'val1', 'key2', 'val2' ...)
% Inputs:
%         data    = [float array] 2-D data array of size (times,trials)
%         srate   = sampling rate
%
% Optional inputs:
%       'cycles'  =  [real positive scalar] Number of cycles in each Morlet
%                   wavelet, constant across frequencies.
%                   or [cycles cycles(2)] wavelet cycles increase with
%                   frequency starting at cycles(1) and,
%                   if cycles(2) > 1, increasing to cycles(2) at
%                   the upper frequency,
%                   or if cycles(2) = 0, same window size at all
%                   frequencies (similar to FFT if cycles(1) = 1)
%                   or if cycles(2) = 1, not increasing (same as giving
%                   only one value for 'cycles'). This corresponds to pure
%                   wavelet with the same number of cycles at each frequencies
%                   if 0 < cycles(2) < 1, linear variation in between pure
%                   wavelets (1) and FFT (0). The exact number of cycles
%                   at the highest frequency is indicated on the command line.
%       'wletmethod' = ['dftfilt2'|'dftfilt3'] Wavelet method/program to use.
%                   {default: 'dftfilt3'}
%                   'dftfilt'  DEPRECATED. Method used in regular timef()
%                              program. Not available any more.
%                   'dftfilt2' Morlet-variant or Hanning DFT (calls dftfilt2()
%                              to generate wavelets).
%                   'dftfilt3' Morlet wavelet or Hanning DFT (exact Tallon
%                              Baudry). Calls dftfilt3().
%                   'comptype'   = ['spectrogram'|'complex'] Compute either the average spectrogram (fast), or the complex wavelet transfor of each frame.
% Optional FFT/DFT parameters:
%      'winsize'  = If cycles==0 (FFT, see 'wavelet' input): data subwindow
%                   length (fastest, 2^n<frames);
%                   if cycles >0: *longest* window length to use. This
%                   determines the lowest output frequency  {~frames/8}
%      'freqs'    = [min max] frequency limits. Default [minfreq srate/2],
%                   minfreq being determined by the number of data points,
%                   cycles and sampling frequency. Enter a single value
%                   to compute spectral decompisition at a single frequency
%                   (note: for FFT the closest frequency will be estimated).
%                   For wavelet, reducing the max frequency reduce
%                   the computation load.
%      'padratio' = FFTlength/winsize (2^k)                     {def: 2}
%                   Multiplies the number of output frequencies by
%                   dividing their spacing. When cycles==0, frequency
%                   spacing is (low_frequency/padratio).
%      'nfreqs'   = number of output frequencies. For FFT, closest computed
%                   frequency will be returned. Overwrite 'padratio' effects
%                   for wavelets. Default: use 'padratio'.
%     'freqscale' = ['log'|'linear'] frequency scale. Default is 'linear'.
%                   Note that for obtaining 'log' spaced freqs using FFT,
%                   closest correspondant frequencies in the 'linear' space
%                   are returned.
%     'wletmethod'= ['dftfilt2'|'dftfilt3'] Wavelet method/program to use.
%                   Default is 'dftfilt3'
%                   'dftfilt3' Morlet wavelet or Hanning DFT
%                   'dftfilt2' Morlet-variant or Hanning DFT.
%                   Note that there are differences betweeen the Hanning
%                   DFTs in the two programs.
%
% Optional time warping:
% Outputs:
%         tf      = complex time frequency array for all trials (freqs,
%                   times, trials)
%         freqs   = vector of computed frequencies (Hz)
%         times   = vector of computed time points (ms)
%         angle   = computed circular mean phase (when evaluated)
%

if nargin < 2
    help timefreq;
    return;
end;

[chan frame trials]= size(data);tmpangle=[];
if trials == 1 && chan ~= 1
    trials = frame;
    frame  = chan;
    chan   = 1;
end;
g = finputcheck(varargin, ...
    {     'winsize'       'integer'  [0 Inf]                  []; ...
    'verbose'       'string'   {'on','off'}             'on'; ...
    'freqs'         'real'     [0 Inf]                  []; ...
    'nfreqs'        'integer'  [0 Inf]                  []; ...
    'freqscale'     'string'   { 'linear','log','' }    'linear'; ...
    'wavelet'       'real'     [0 Inf]                  0; ...
    'cycles'        {'real','integer'}    [0 Inf]       4; ...
    'padratio'      'integer'  [1 Inf]                  2; ...
    'subitc'        'string'   {'on','off'}             'off'; ...
    'wletmethod'    'string'   {'dftfilt2','dftfilt3'}    'dftfilt3'; ...
    'comptype' 'string'   {'spectrogram','complex'}         'complex'; ...
    });
if isstr(g), error(g); end;
if isempty(g.freqscale), g.freqscale = 'linear'; end;
if isempty(g.winsize),   g.winsize   = max(pow2(nextpow2(frame)-3),4); end;
if isempty(g.freqs),     g.freqs     = [0 srate/2]; end;

% checkin parameters
% ------------------

% Use 'wavelet' if 'cycles' undefined for backwards compatibility
if g.cycles == 0
    g.cycles = g.wavelet;
end

if (g.winsize > frame)
    error('Value of winsize must be less than frame length.');
end
if (pow2(nextpow2(g.padratio)) ~= g.padratio)
    error('Value of padratio must be an integer power of two [1,2,4,8,16,...]');
end

% finding frequency limits
% ------------------------
if g.cycles(1) ~= 0 & g.freqs(1) == 0, g.freqs(1) = srate*g.cycles(1)/g.winsize; end;

% finding frequencies
% -------------------
if length(g.freqs) == 2
    
    % min and max
    % -----------
    if g.freqs(1) == 0 & g.cycles(1) ~= 0
        g.freqs(1) = srate*g.cycles(1)/g.winsize;
    end;
    
    % default number of freqs using padratio
    % --------------------------------------
    if isempty(g.nfreqs)
        g.nfreqs = g.winsize/2*g.padratio+1;
        % adjust nfreqs depending on frequency range
        tmpfreqs = linspace(0, srate/2, g.nfreqs);
        tmpfreqs = tmpfreqs(2:end);  % remove DC (match the output of PSD)
        
        
        
        % find number of frequencies
        % --------------------------
        g.nfreqs = length(tmpfreqs( intersect( find(tmpfreqs >= g.freqs(1)), find(tmpfreqs <= g.freqs(2)))));
        if g.freqs(1)==g.freqs(2), g.nfreqs = 1; end;
    end;
    
    % find closest freqs for FFT
    % --------------------------
    if strcmpi(g.freqscale, 'log')
        g.freqs = linspace(log(g.freqs(1)), log(g.freqs(end)), g.nfreqs);
        g.freqs = exp(g.freqs);
    else
        g.freqs = linspace(g.freqs(1), g.freqs(2), g.nfreqs); % this should be OK for FFT
        % because of the limit adjustment
    end;
end;
g.nfreqs = length(g.freqs);

% function for time freq initialisation
% -------------------------------------

freqs = g.freqs;
if length(g.cycles) == 2
    if g.cycles(2) < 1
        g.cycles = [ g.cycles(1) g.cycles(1)*g.freqs(end)/g.freqs(1)*(1-g.cycles(2))];
    end
    verboseprintf(g.verbose, 'Using %g cycles at lowest frequency to %g at highest.\n', g.cycles(1), g.cycles(2));
elseif length(g.cycles) == 1
    verboseprintf(g.verbose, 'Using %d cycles at all frequencies.\n',g.cycles);
else
    verboseprintf(g.verbose, 'Using user-defined cycle for each frequency\n');
end
if strcmp(g.wletmethod, 'dftfilt2')
    g.win    = dftfilt2(g.freqs,g.cycles,srate, g.freqscale); % uses Morlet taper by default
elseif strcmp(g.wletmethod, 'dftfilt3')     % Default
    g.win    = dftfilt3(g.freqs,g.cycles,srate, 'cycleinc', g.freqscale); % uses Morlet taper by default
else return
end
g.winsize = 0;
for index = 1:length(g.win)
    g.winsize = max(g.winsize,length(g.win{index}));
end;

% -------------------------------
% compute time freq decomposition
% -------------------------------
verboseprintf(g.verbose, 'The window size used is %d samples (%g ms) wide.\n',g.winsize, 1000/srate*g.winsize);

% % prepare wavelet filters
% % -----------------------
% for index = 1:length(g.win)
%     g.win{index} = transpose(repmat(g.win{index}, [trials 1]));
% end;

% apply filters
% -------------

maxpad=cellfun(@(x) length(x),g.win);
if 0
    %modify order of data dimensions: time first  %MB local change
    data=permute(data,[2 3 1]);
end
datapad=cat(1,zeros(max(maxpad),size(data,2),size(data,3)),data,zeros(max(maxpad),size(data,2),size(data,3)));
datapad=permute(datapad,[1 3 2]);
datafft=fft(datapad);

verboseprintf(g.verbose,'.');

sigma=g.cycles;

%allow for varying number of cycles
scale=1./(2*pi*freqs./sigma);

%evalute analytic experession of the Fourier transform of the wavelet
csigma=1/sqrt(1+exp(-sigma^2)-2*exp(-3/4*sigma.^2));
ksigma=exp(-1/2*sigma^2);
freqind=(0:(size(datapad,1)-1))/size(datapad,1)*srate;
wavtf=csigma*pi^(-1/4)*sqrt(ones(length(freqind),1)*scale).*(exp(-1/2*(sigma-2*pi*freqind'*scale).^2)-ksigma*exp(-1/2*(2*pi*freqind'*scale).^2));

switch g.comptype
    case 'complex'
        %convres=ifft((repmat(datafft,[1,length(freqs),1])).*(repmat(wavtf,[1,1,size(datafft,3)])));
        convres=ifft(bsxfun(@times,datafft,wavtf));
        
        tmpall = permute(convres((max(maxpad)+1):(max(maxpad)+size(data,1)),:,:),[2 1 3]);
    case 'spectrogram'% faster! compute directly average power and circular mean across trials
        convres=0;convangle=0;
        for ktrial=1:size(datafft,3)
            %tmpconv=ifft((repmat(datafft(:,:,ktrial),[1,length(freqs),1])).*wavtf);
            tmpconv=ifft(bsxfun(@times,datafft(:,:,ktrial),wavtf));
            convres=convres+abs(tmpconv).^2;
            convangle=convangle+tmpconv./abs(tmpconv);
        end
        convres=convres/size(datafft,3);
        convangle=convangle/size(datafft,3);
        tmpall = permute(convres((max(maxpad)+1):(max(maxpad)+size(data,1)),:,:),[2 1 3]);
        tmpangle=permute(convangle((max(maxpad)+1):(max(maxpad)+size(data,1)),:,:),[2 1 3]);
end

timesout=(0:(size(data,1)-1))/srate*1000;
verboseprintf(g.verbose, '\n');
return;



function verboseprintf(verbose, varargin)
if strcmpi(verbose, 'on')
    fprintf(varargin{:});
end;

