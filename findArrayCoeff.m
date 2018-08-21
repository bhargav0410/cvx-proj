function [w_1] = findArrayCoeff(alp, bet, num_ants, phase_rot, user_doa, int_doa, beamwidth, corrMat, array_of_angles, type)
    
for p = array_of_angles
    v_vec(p,:) = [real(phase_rot(p,:)),imag(phase_rot(p,:))];
    v_cap(p,:) = [-imag(phase_rot(p,:)),real(phase_rot(p,:))];
    v_mat(p,:,:) = [real(phase_rot(p,:)),imag(phase_rot(p,:));-imag(phase_rot(p,:)),real(phase_rot(p,:))];
end 
    corr_mat = corrMat;
    % MVDR only
    if type == 1
        cvx_begin quiet
            variables w(num_ants*2)
            temp(:,:) = sqrtm(corr_mat(:,:));
            temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
            min_val = [temp_*w];
            minimize norm(min_val,2)
                subject to
                    v_vec(user_doa,:)*w == num_ants;
                    v_cap(user_doa,:)*w == 0;
                       
        cvx_end
        w_1 = w(1:num_ants) + 1i*w(num_ants+1:end);
    % MVDR with sidelobe suppression
    elseif type == 2
        cvx_begin quiet
            variables w(num_ants*2) e
            temp(:,:) = sqrtm(corr_mat(:,:));
            temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
            min_val = [temp_*w;alp*e];
            minimize norm(min_val,2)
                subject to
                    for p = array_of_angles
                        if p == user_doa
                            v_vec(p,:)*w == num_ants;
                            v_cap(p,:)*w == 0;
                        elseif norm(p - user_doa,1) >= beamwidth/2
                            temp_v(:,:) = v_mat(p,:,:);
                            norm(temp_v*w,2) <= e;
                        end
                    end
        cvx_end
        w_1 = w(1:num_ants) + 1i*w(num_ants+1:end);
    
    % MVDR with sidelobe and interference suppression
    elseif type == 3
        alp = 1;
        cvx_begin quiet
            variables w(num_ants*2) e u
            temp(:,:) = sqrtm(corr_mat(:,:));
            temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
            min_val = [temp_*w;alp*e;bet*u];
            minimize norm(min_val,2)
                subject to
                    for p = array_of_angles
                        if p == user_doa
                            v_vec(p,:)*w == num_ants;
                            v_cap(p,:)*w == 0;
                        elseif ismember(p,int_doa)
                            temp_v(:,:) = v_mat(p,:,:);
                            norm(temp_v*w,2) <= u;
                        elseif norm(p - user_doa,1) >= beamwidth/2
                            temp_v(:,:) = v_mat(p,:,:);
                            norm(temp_v*w,2) <= e;
                        end
                    end
        cvx_end
        w_1 = w(1:num_ants) + 1i*w(num_ants+1:end);
    
    % MVDR with interference suppression
    elseif type == 4
        alp = 1;
        cvx_begin quiet
            variables w(num_ants*2) u
            temp(:,:) = sqrtm(corr_mat(:,:));
            temp_ = [real(temp),-imag(temp);imag(temp),real(temp)];
            min_val = [temp_*w;bet*u];
            minimize norm(min_val,2)
                subject to
                    for p = array_of_angles
                        if p == user_doa
                            v_vec(p,:)*w == num_ants;
                            v_cap(p,:)*w == 0;
                        elseif ismember(p,int_doa)
                            temp_v(:,:) = v_mat(p,:,:);
                            norm(temp_v*w,2) <= u;
                        end
                    end
        cvx_end
        w_1 = w(1:num_ants) + 1i*w(num_ants+1:end);
    end
end