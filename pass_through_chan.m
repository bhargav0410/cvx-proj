function [y_rcr] = pass_through_chan(Fs_upsampled, rate, Y, numTx, numRx, path_delays, path_gains, snr, filters)

y_rcr = [];

% Raised Cosine Transmit Filter

if filters == 1
    Nsym = 6;           % Filter span in symbol durations
    beta = 0.1;         % Roll-off factor
    sampsPerSym = Fs_upsampled/rate;    % Upsampling factor

    rctFilt = comm.RaisedCosineTransmitFilter(...
      'Shape',                  'Normal', ...
      'RolloffFactor',          beta, ...
      'FilterSpanInSymbols',    Nsym, ...
      'OutputSamplesPerSymbol', sampsPerSym);

    y_rct = rctFilt(Y);
else
    y_rct = Y;
end
%%%% If MIMO channel is to be used %%%%
h = comm.MIMOChannel;
h.SampleRate = Fs_upsampled;
h.SpatialCorrelation = false; % Independent channels
h.NumTransmitAntennas = numTx;
h.NumReceiveAntennas = numRx;
h.FadingDistribution = 'Rayleigh';
% h.DirectPathDopplerShift = Doppshift;
% h.MaximumDopplerShift = Doppshift;
% h.KFactor = 5;
h.PathDelays = path_delays; 
h.NormalizePathGains = false;
h.AveragePathGains = path_gains;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_filt = step(h,y_rct);

% Raised Cosine Receive Filter

if filters == 1
    Nsym = 6;           % Filter span in symbol durations
    beta = 0.1;         % Roll-off factor
    sampsPerSym = Fs_upsampled/rate;    % Downsampling factor

    rcrFilt = comm.RaisedCosineReceiveFilter(...
      'Shape',                  'Normal', ...
      'RolloffFactor',          beta, ...
      'FilterSpanInSymbols',    Nsym, ...
      'InputSamplesPerSymbol', sampsPerSym, ...
      'DecimationFactor', sampsPerSym );

    y_rcr = rcrFilt(y_filt);
else
    y_rcr = y_filt;
end

% % Adding noise
% y_rcr = awgn(y_rcr,snr);

end