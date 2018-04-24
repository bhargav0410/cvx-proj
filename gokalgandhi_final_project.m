clc;
clear all;
close all;

Fc = 5400e6;
Fs = 10e6;
num_ants = 16;
lmd = 3e8/Fc;
d = lmd/2;
k_ = 2*pi/lmd;
theta = 90;
phi = 1:360;
e = 0.001;

user_doa = 270;

% OFDM packet generation
num_syms = 1000;
k = 1;
NFFT = 64;
cp_len = 16;

% preamble = getOFDMPreambleAndPilot('Preamble', 64, [6; 5]);
% preamble = preamble/max(abs(preamble));
% ofdmMod = comm.OFDMModulator('FFTLength',NFFT,'CyclicPrefixLength',cp_len,'NumGuardBandCarriers',[0;0],'InsertDCNull',true,'NumSymbols',num_syms);
% ofdm_demod = comm.OFDMDemodulator('FFTLength',NFFT,'CyclicPrefixLength',cp_len,'NumGuardBandCarriers',[0;0],'RemoveDCCarrier',true,'NumSymbols',num_syms);
% X_Input = randi([0 2^k - 1],(NFFT-1)*ofdmMod.NumSymbols,1);
% YQPSKOutput = qammod(X_Input,2^k,'UnitAveragePower',true);
% Y_QPSKOutput = reshape(YQPSKOutput,NFFT-1,ofdmMod.NumSymbols);
% Y_OFDMOutput = step(ofdmMod, Y_QPSKOutput);%/(max(real(YQPSKOutput)/9)));
% OFDMOutput = [reshape(Y_OFDMOutput,(NFFT+cp_len)*ofdmMod.NumSymbols,1)/max(abs(Y_OFDMOutput))];
% 
% X_Input_int = randi([0 2^k - 1],(NFFT-1)*ofdmMod.NumSymbols,1);
% YQPSKOutput_int = qammod(X_Input,2^k,'UnitAveragePower',true);
% Y_QPSKOutput_int = reshape(YQPSKOutput,NFFT-1,ofdmMod.NumSymbols);
% Y_OFDMOutput_int = step(ofdmMod, Y_QPSKOutput);%/(max(real(YQPSKOutput)/9)));
% OFDMOutput_int = [reshape(Y_OFDMOutput,(NFFT+cp_len)*ofdmMod.NumSymbols,1)/max(abs(Y_OFDMOutput))];
% 
% % MIMO channel
% mimo_chan = comm.MIMOChannel;
% mimo_chan.SampleRate = Fs;
% mimo_chan.PathDelays = [0 1 2 3 4 5]*1e-7;
% mimo_chan.AveragePathGains = [0 -10 -20 -30 -40 -50];
% mimo_chan.FadingDistribution = 'Rician';
% mimo_chan.SpatialCorrelationSpecification = 'None';
% mimo_chan.NumTransmitAntennas = 1;
% mimo_chan.NumReceiveAntennas = num_ants;
% 
% % ULA example
% for p = phi
%     for m  = 1:num_ants
%         phase_rot(p,m) = exp(1i * (m-1)*k_*d*cos(p*pi/180));
%     end
% end
% 
% % UPA example
% % for p = phi
% %     for m  = 1:sqrt(num_ants)
% %         for n = 1:sqrt(num_ants)
% %             phase_rot(p,(m-1)*sqrt(num_ants) + n) = exp(1i * (m-1)*k_*d*cos(p*pi/180))*exp(1i * (n-1)*k_*d*cos(p*pi/180));
% %         end
% %     end
% % end
% % for m  = 1:num_ants
% %     phase_rot(m) = exp(1i * (m-1)*k*d*sin(theta)*cos(phi));
% % end
% 
% % Reception at antenna array and addition of noise
% rec_packet = mimo_chan(OFDMOutput).*phase_rot(user_doa,:) + OFDMOutput_int.*phase_rot(310,:);% + OFDMOutput_int.*phase_rot(20,:) + OFDMOutput_int.*phase_rot(330,:);
% 
% % for m = 1:num_ants
% %     Y(:,:,m) = ofdm_demod(awgn(rec_packet(:,m),3,'measured'));
% % end
% for m = 1:num_ants
%     Y(:,m) = awgn(rec_packet(:,m),20,'measured');
% end

% for i = 1:NFFT-1
%     temp_y(:,:) = Y(i,1:min(2,num_syms),:);
%     corr_mat(:,:,i) = 1/(min(2,num_syms)) * transpose(temp_y)*conj(temp_y);
% end

% For use durring testbed experimentation only !!!!!

for i = 1:num_ants
    rec(:,i) = read_float_binary('D:\Massive_MIMO_Programs\corr_rec_ch_',char(num2str(i-1)),'_binary');
end


corr_mat = (1/(NFFT-1)) * (Y(1:NFFT-1,:)'*Y(1:NFFT-1,:));

for p = phi
    v_vec(p,:) = ([real(phase_rot(p,:)),imag(phase_rot(p,:))]);
    v_cap(p,:) = ([-imag(phase_rot(p,:)),real(phase_rot(p,:))]);
    v_mat(p,:,:) = [real(phase_rot(p,:)),imag(phase_rot(p,:));-imag(phase_rot(p,:)),real(phase_rot(p,:))];
end
    
% cvx_begin quiet
%     variables w_max(num_ants*2)
%         temp(:,:) = sqrtm(corr_mat(:,:));
%         temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
%         min_val = [temp_*w_max];
%         minimize norm(min_val,2)
% %         subject to
% %                 norm(temp_*w) <= t
% %             v_vec(user_doa,:)*w <= num_ants
% %             v_vec(user_doa,:)*w == 1
% %             v_cap(user_doa,:)*w == 0
% %             for p = [180:250,290:360]
% %                 temp_v(:,:) = v_mat(p,:,:);
% %                 norm(temp_v*w,2) <= e
% %             end
% cvx_end

% [eig_vec,eig_val] = eig(corr_mat);
% for i = 1:360
%     P(i) = 1/(conj(phase_rot(i,:))*eig_vec(:,3:num_ants)*eig_vec(:,3:num_ants)'*transpose(phase_rot(i,:)));
% end
% figure;plot(180:359,abs(P(180:359)));
% [max_val,max_arg] = max(abs(P(180:359)));
% user_doa = max_arg+179;

% cvx_solver sedumi
alp = 1e2;i=1;
figure;
for alp = 1e6
    cvx_begin quiet
        variables w(num_ants*2) e
            temp(:,:) = sqrtm(corr_mat(:,:,i));
            temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
            min_val = [temp_*w;alp*e];
            minimize norm(min_val,2)
            subject to
%                 norm(temp_*w) <= t
                v_vec(user_doa,:)*w <= num_ants
                v_vec(user_doa,:)*w >= 1
                v_cap(user_doa,:)*w == 0
                for p = [180:user_doa-20,user_doa+20:360]
                    temp_v(:,:) = v_mat(p,:,:);
                    norm(temp_v*w,2) <= e
                end
    cvx_end

%     cvx_begin quiet
%         variables w(num_ants*2) t
%         temp(:,:) = sqrtm(corr_mat(:,:,i));
%         temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
%         minimize t
%         subject to
%             norm(temp_*w,2) <= t
%             v_vec(user_doa,:)*w <= num_ants
%             1 - v_vec(user_doa,:)*w <= 0
%             v_cap(user_doa,:)*w == 0
%             for p = [25:90,270:345]
%                 temp_v(:,:) = v_mat(p,:,:);
%                 norm(temp_v*w,2) <= sqrt(e)
%             end
%     cvx_end
            
        
    w_(:,i) = w(1:num_ants) + 1i*w(num_ants+1:end);
    plot(1:360,(20*log10(abs(phase_rot([1:360],:)*conj(w_))))); hold on;
end


% scatterplot(Y*conj(w_));
for i = 1:num_ants
    demod_op(:,:,i) = ofdm_demod(Y(:,i)*conj(w_(i)));
end
rec = sum(demod_op,3);
chan_est = rec(:,1)./Y_QPSKOutput(:,1);
scatterplot(reshape(rec(:,2:end)./repmat(chan_est,1,num_syms-1),(NFFT-1)*(num_syms-1),1));

QPSK_out = reshape(rec(:,2:end)./repmat(chan_est,1,num_syms-1),(NFFT-1)*(num_syms-1),1);
% for i = 1:num_syms
%     temp_demod(:,:) = Y(:,i,:);
%     op(:,i) = sum(temp_demod.*w_',2);
QPSK_demod = qamdemod(QPSK_out*2./max(abs(QPSK_out)),2^k,'UnitAveragePower',true);

err = sum(QPSK_demod ~= X_Input(NFFT:end))/length(X_Input(NFFT:end));
% end
% chan_est = op(:,1)./Y_QPSKOutput(:,1);
% rec = op(:,2:end)./repmat(chan_est,1,num_syms-1);
% scatterplot(reshape(op,(NFFT-1)*num_syms,1));