function [w,t] = powerOptimal(H, eig_val, eig_vec, N_bar, gamma)
    M = 16; %no. of users
    G = 3; %no. of user groups
    K = G;
    %g{1} = [1,4,6,9,16]; g{2} = [2,7,8,11,14]; g{3} = [3,5,10,12,13,15];
    g{1} = [1,2,3,4,5]; g{2} = [6,7,8,9,10]; g{3} = [11,12,13,14,15,16];
    N = 6; %no. of transmit antennas
    %gamma = .001; %QoS

    %for a = 1:length(ang_spread)


    %%power minimization
    for ml = 1:4

        cvx_begin quiet

        variables obj;
        variables w(G*N+1,1) Complex;

        expression s(M,1)
        expression sr(M,1)
        for k = 1:K
            for ma = 1:length(g{k})
                m = g{k}(ma);
                t = 0;
                for n = 1:length(N_bar(m))
                    eval = eig_val{m};
                    evec = eig_vec{m};
                    t = t + sqrt(eval(n))*real(w((N*(k-1))+2:(N*(k-1))+1+N)'*evec(:,n));
                end
                s(m) = (1/sqrt(N_bar(m)))*t;
                sr(m) = sqrt(gamma)*norm(H{m,k}*w,2);%sum_square_abs(H{m,k}*w);
            end
        end

        minimize obj;
        subject to
        obj >= norm(w,2);%sum_square_abs(w);
        w(1) == 1;
        for m = 1:M
            s(m) >= sr(m);
        end
        cvx_end
        clear s sr;

        %% descending eigen values & eigen vectors
        dec_vec = cell(M,1);
        for k = 1:K
            for ma = 1:length(g{k})
                m = g{k}(ma);
                for n = 1:length(eig_vec{m})
                    dec_vec{m}(n) = eig_val{m}(n)*abs(w((N*(k-1))+2:(N*(k-1))+1+N)'*eig_vec{m}(:,n))^2;
                end
            end
        end
        for m = 1:M
            [Y,I] = sort(dec_vec{m},'descend');
            eig_val{m} = eig_val{m}(I);
            eig_vec{m} = eig_vec{m}(I,:);
        end

        %% finding optimum set%
        for k = 1:K
            for ma = 1:length(g{k})
                m = g{k}(ma);
                for n = length(eig_val{m}):-1:1
                    ab = 0;
                    for nn = 1:n
                        ab = ab + sqrt(eig_val{m}(nn))*abs(w((N*(k-1))+2:(N*(k-1))+1+N)'*eig_vec{m}(:,nn));
                    end
                    arg_val{m}(n) = 1/sqrt(n)*ab;
                end
                [Mag,new_n] = max(arg_val{m});
                N_bar(m) = new_n;
            end
        end

        %% rotation approximation%
        for k = 1:K
            for ma = 1:length(g{k})
                m = g{k}(ma);
                for n = 1:length(eig_vec{m})
                    alpha = angle(w((N*(k-1))+2:(N*(k-1))+1+N)'*eig_vec{m}(:,n));
                    eig_vec{m}(:,n) = eig_vec{m}(:,n)*exp(-1i*alpha);
                end
            end
        end
    end
end
