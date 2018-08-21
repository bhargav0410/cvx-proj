clc;
clear all;
close all;

% This code generates OFDM symbols for intended and interference source of
% signals received on an array, then AWGN is added. Then the correlation
% matrix is calculated for the received signal, and then optimum weights
% for MVDR, MVDR with sidelobe suppression and MVDR with sidelobe and
% interference suppression are calculated using the SeDuMi toolbox. The
% resulting beampatterns are then plotted for all three cases.

% Simulation variables
Fc = 5400e6; % Center frequency
Fs = 10e6; % Sampling rate
num_ants = 16; % Number of antennas in the array
lmd = 3e8/Fc; % Wavelength
d = lmd/2; % Distance between antennas (for simulation)
k_ = 2*pi/lmd; % Wavenumber
theta = 90; % Elevation angle
phi = 1:360; % Azimuth angle 

% Angle of arrival for intended user in simulation
user_doa = 90;

% OFDM packet generation
k = 1; % Number of bits per symbol
NFFT = 64; % FFT size in OFDM
cp_len = 16; % Cyclic prefix length
num_syms = ceil((1e5)/(NFFT-1)); % Number of symbols
% 
% Defining the OFDM modulator and demodulator
ofdmMod = comm.OFDMModulator('FFTLength',NFFT,'CyclicPrefixLength',cp_len,'NumGuardBandCarriers',[0;0],'InsertDCNull',true,'NumSymbols',num_syms);
ofdm_demod = comm.OFDMDemodulator('FFTLength',NFFT,'CyclicPrefixLength',cp_len,'NumGuardBandCarriers',[0;0],'RemoveDCCarrier',true,'NumSymbols',num_syms);

% % Formation of OFDM symbol for intended source
% X_Input = randi([0 2^k - 1],(NFFT-1)*ofdmMod.NumSymbols,1); % Generation of random symbols
% YBPSKOutput = qammod(X_Input,2^k,'UnitAveragePower',true); % Modulation into BPSK
% Y_BPSKOutput = reshape(YBPSKOutput,NFFT-1,ofdmMod.NumSymbols); 
% Y_OFDMOutput = step(ofdmMod, Y_BPSKOutput); % Conversion into an OFDM frame
% OFDMOutput = [reshape(Y_OFDMOutput,(NFFT+cp_len)*ofdmMod.NumSymbols,1)/max(abs(Y_OFDMOutput))];
% 
% % Formation of OFDM symbol for interference source 
% X_Input_int = randi([0 2^k - 1],(NFFT-1)*ofdmMod.NumSymbols,1);
% YBPSKOutput_int = qammod(X_Input_int,2^k,'UnitAveragePower',true);
% Y_BPSKOutput_int = reshape(YBPSKOutput_int,NFFT-1,ofdmMod.NumSymbols);
% Y_OFDMOutput_int = step(ofdmMod, Y_BPSKOutput_int);
% OFDMOutput_int = [reshape(Y_OFDMOutput_int,(NFFT+cp_len)*ofdmMod.NumSymbols,1)/max(abs(Y_OFDMOutput))];

% Calculation of steering vectors for ULA example
for p = phi
    for m  = 1:num_ants
        phase_rot(p,m) = exp(1i * (m-1)*k_*d*cos(p*pi/180));
    end
end

% Calculation of steering vecotrs for UPA example (not shown in report)
% for p = phi
%     for m  = 1:sqrt(num_ants)
%         for n = 1:sqrt(num_ants)
%             phase_rot(p,(m-1)*sqrt(num_ants) + n) = exp(1i * (m-1)*k_*d*cos(p*pi/180))*exp(1i * (n-1)*k_*d*sin(p*pi/180));
%         end
%     end
% end

% Testbed array
xpos = [0 0.0508*2 0.0508*4 0.0508*6 0.0508 0.0508*3 0.0508*5 0.0508*7 0 0.0508*2 0.0508*4 0.0508*6 0.0508 0.0508*3 0.0508*5 0.0508*7];
ypos = [0 0 0 0 0.0508 0.0508 0.0508 0.0508 0.0508*2 0.0508*2 0.0508*2 0.0508*2 0.0508*3 0.0508*3 0.0508*3 0.0508*3];
figure;
plot(xpos,ypos, '*');
axis([-0.05 0.4 -0.05 0.2]);
xlabel('Distance (m)');
ylabel('Distance (m)');

% Calculation of steering vectors for the testbed arrays
for i = phi
    for m = 1:num_ants
            phase_rot(i,m) = exp(1i*2*pi*(1/lmd)*(cos(i*pi/180)*xpos(m) + sin(i*pi/180)*ypos(m)));
    end
end


% Reception at antenna array and addition of noise
% rec_packet = OFDMOutput.*phase_rot(user_doa,:) + OFDMOutput_int.*phase_rot(60,:);
% for m = 1:num_ants
%     Y(:,m) = reshape(ofdm_demod(awgn(rec_packet(:,m),3,'measured')),(NFFT-1)*num_syms,1);
% end


% --------------- For use during testbed experimentation only !!!!!-------

for i = 1:num_ants
    rec = read_float_binary(['C:\Users\Bhargav04\Documents\Massive MIMO programs\corr_rec_2_iter2_ch_',char(num2str(i-1)),'_binary']);
    Y(:,i) = rec(1:2:end) + 1i*rec(2:2:end);
end
rec = read_float_binary(['C:\Users\Bhargav04\Documents\Massive MIMO programs\OFDM_BPSK_64FFT\OFDMPacket.dat']);
X = rec(1:2:end) + 1i*rec(2:2:end);
X_QPSK = ofdm_demod(X(256:end));
for i = 1:num_ants
    rec = read_float_binary(['C:\Users\Bhargav04\Documents\Massive MIMO programs\noise_ch_',char(num2str(i-1)),'_binary']);
    noise(:,i) = rec(1:2:end) + 1i*rec(2:2:end);
end
noise_mean = max(mean(abs(noise).^2));
% 
% ------------------------------------------------------------------------


% Calculation of correlation matrix
corr_iter = 10000;
corr_mat = (1/(corr_iter)) * (Y(1:corr_iter,:)'*Y(1:corr_iter,:));


% --------------- For use during testbed experimentation only !!!!!-------
[eig_vec,eig_val] = eig(corr_mat);
for i = 1:360
%     P(i) = 1/(conj(phase_rot(i,:))*eig_vec(:,3:num_ants)*eig_vec(:,3:num_ants)'*transpose(phase_rot(i,:)));
    P(i) = 1/(conj(phase_rot(i,:))*eig_vec(:,max(2,sum(diag(eig_val)>noise_mean)+1):num_ants)*eig_vec(:,max(2,sum(diag(eig_val)>noise_mean)+1):num_ants)'*transpose(phase_rot(i,:)));
end
figure;plot(1:360,abs(P));
[max_val,max_arg] = max(abs(P));
user_doa = max_arg;
% ------------------------------------------------------------------------


% Divding complex numbers into real and imaginary so that the toolbox can
% use the data
for p = phi
    v_vec(p,:) = [real(phase_rot(p,:)),imag(phase_rot(p,:))];
    v_cap(p,:) = [-imag(phase_rot(p,:)),real(phase_rot(p,:))];
    v_mat(p,:,:) = [real(phase_rot(p,:)),imag(phase_rot(p,:));-imag(phase_rot(p,:)),real(phase_rot(p,:))];
end
temp(:,:) = sqrtm(corr_mat(:,:));
temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];

% Using SeDuMi for finding the optimized weights
cvx_solver sedumi
% Using very high values for penalty function weights
alp = 1e6;
bet = 1e6;
i = 1;

figure;

% Use the for loop below to change the values of alpha for variation in the
% sidelobe threshold. In the report values from 1 to 50 were used
%{
for alp = 1e6
    % MVDR only
    cvx_begin quiet
        variables w(num_ants*2)
            min_val = [temp_*w];
            min_val2 = alp*[];
            minimize norm(min_val,2)
            subject to
                v_vec(user_doa,:)*w == num_ants
                v_cap(user_doa,:)*w == 0
    cvx_end
    w_1(:,i) = w(1:num_ants) + 1i*w(num_ants+1:end);
    e_1 = max(abs(phase_rot([1:user_doa-20,user_doa+20:360],:)*conj(w_1)));
    
    % MVDR with sidelobe suppression
    cvx_begin %quiet
        variables w(num_ants*2) e
%             temp(:,:) = sqrtm(corr_mat(:,:,i));
%             temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
            min_val = [temp_*w;alp*e];
            minimize norm(min_val,2)
            subject to
                v_vec(user_doa,:)*w == num_ants
                v_cap(user_doa,:)*w == 0
                for p = [1:user_doa-20,user_doa+20:180]
                    temp_v(:,:) = v_mat(p,:,:);
                    norm(temp_v*w,2) <= e
                end
    cvx_end
    e_2 = e;
    w_2(:,i) = w(1:num_ants) + 1i*w(num_ants+1:end);
    
    % MVDR with sidelobe and interference suppression
    alp = 1;
    cvx_begin %quiet
        variables w(num_ants*2) e u
%             temp(:,:) = sqrtm(corr_mat(:,:,i));
%             temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
            min_val = [temp_*w;alp*e;bet*u];
            minimize norm(min_val,2)
            subject to
                v_vec(user_doa,:)*w == num_ants
                v_cap(user_doa,:)*w == 0
                v_vec(126,:)*w == num_ants
                v_cap(126,:)*w == 0
                v_vec(239,:)*w == num_ants
                v_cap(239,:)*w == 0
                v_vec(307,:)*w == num_ants
                v_cap(307,:)*w == 0
%                 for p = [1:user_doa-20,user_doa+20:180]
%                     temp_v(:,:) = v_mat(p,:,:);
%                     norm(temp_v*w,2) <= e
%                 end
                for p = [240]
                    temp_v(:,:) = v_mat(p,:,:);
                    norm(temp_v*w,2) <= u;
                end
    cvx_end
    e_3 = e;
    w_3(:,i) = w(1:num_ants) + 1i*w(num_ants+1:end);
    
    % PLotting for all three cases above
    plot((1:360),(10*log10(abs(phase_rot(1:360,:)*conj(w_1)))),'r'); hold on;
    plot((1:360),(10*log10(abs(phase_rot(1:360,:)*conj(w_2)))),'b--'); hold on;
    plot((1:360),(10*log10(abs(phase_rot(1:360,:)*conj(w_3)))),'k:'); hold on;
%     eps_(alp) = e; % To be uncommented when plotting the variation in
%     sidelobe levels w.r.t alpha
end
%}

w(1,:) = findArrayCoeff(alp, bet, num_ants, phase_rot, user_doa, 239, 30, corr_mat, phi, 1);
w(2,:) = findArrayCoeff(alp, bet, num_ants, phase_rot, user_doa, 239, 30, corr_mat, phi, 2);
w(3,:) = findArrayCoeff(alp, bet, num_ants, phase_rot, user_doa, 239, 30, corr_mat, phi, 3);
w(4,:) = findArrayCoeff(alp, bet, num_ants, phase_rot, user_doa, 239, 30, corr_mat, phi, 4);

plot(phi,(10*log10(abs(phase_rot(phi,:)*w(1,:)'))),'r'); hold on;
plot(phi,(10*log10(abs(phase_rot(phi,:)*w(2,:)'))),'b--'); hold on;
plot(phi,(10*log10(abs(phase_rot(phi,:)*w(3,:)'))),'k:'); hold on;
plot(phi,(10*log10(abs(phase_rot(phi,:)*w(4,:)'))),'g-'); hold on;
legend('MVDR only','with sidelobe suppression','with sidelobe and interference suppression','with interference suppression');

grid on;
xlabel('Direction of signal');
ylabel('Beampattern (dB)');
legend('Without sidelobe suppression','With sidelobe suppression','With sidelobe and interference suppression');
title('Beampatterns of the array after optimum weights are calculated for above three cases');
% axis([1 180 -100 0]);

% Plot for variation in alpha
% figure; plot(eps_, 'r-');
% xlabel('\alpha');
% ylabel('Sidelobe level threshold (in dB)');
% title('Variation in sidelobe level threshold \epsilon w.r.t the weight \alpha');



% ------------------ For deocding of OFDM symbols ------------------------
for i = 1:num_ants
    demod_op(:,i) = Y(:,i)*conj(w_3(i));
end

rec = ofdm_demod(sum(demod_op,2));
chan_est = rec(:,1)./X_QPSK(:,1); % Use for LS channel estimation for testbed

% % LS estimation for OFDM frame
% chan_est = rec(:,1)./Y_BPSKOutput(:,1);
% 
QPSK_out = reshape(rec(:,2:end)./repmat(chan_est,1,num_syms-1),(NFFT-1)*(num_syms-1),1);
scatterplot(QPSK_out*2/max(abs(QPSK_out)));
QPSK_demod = qamdemod(QPSK_out*2./max(abs(QPSK_out)),2^k,'UnitAveragePower',true);

% % Use below error calculation for testbed experiments
err = sum(QPSK_demod ~= qamdemod(reshape(X_QPSK(:,2:end)./max(max(abs(X_QPSK(:,2:end)))),(NFFT-1)*(num_syms-1),1),2^k,'UnitAveragePower',true))/length(X_QPSK(NFFT:end));

% % Use below error calculation for BER output
% err = sum(QPSK_demod ~= X_Input(NFFT:end))/length(X_Input(NFFT:end));