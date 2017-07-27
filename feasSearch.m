function [ cvx_optval, z, eig_val, eig_vec ] = feasSearch( g, M, H, eig_val, eig_vec, N_bar, gamma )
G = 3; %no. of user groups
K = G;
N = 6;
for ml = 1:3
        cvx_begin quiet
        variables z(M,1);
        variables obj
        variables w(G*N+1,1) Complex;
        expression s(M,1)
        expression sr(M,1)
        for k = 1:K
            for ma = 1:length(g{k})
                m = g{k}(ma);
                t = 0;
                for n = 1:N_bar(m)
                    eval = eig_val{m};
                    evec = eig_vec{m};
                    t = t + sqrt(eval(n))*real(w((N*(k-1))+2:(N*(k-1))+1+N)'*evec(:,n));
                end
                s(m) = (1/sqrt(N_bar(m)))*t;
                sr(m) = sqrt(gamma)*norm(H{m,k}*w,2);%sum_square_abs(H{m,k}*w);
            end
        end 
        
        minimize ones(1,M)*z;
        subject to
        for m = 1:M
            sr(m) - s(m) <= z(m);
            z(m) >= 0;
        end
        w(1) == 1;
        cvx_end


     %% power minimization
    %for ml = 1:4

    %     cvx_begin quiet
    % 
    %     variables obj;
    %     variables w(G*N+1,1) Complex;
    % 
    %     expression s(M,1)
    %     expression sr(M,1)
    %     for k = 1:K
    %         for ma = 1:length(g{k})
    %             m = g{k}(ma);
    %             t = 0;
    %             for n = 1:length(N_bar(m))
    %                 eval = eig_val{m};
    %                 evec = eig_vec{m};
    %                 t = t + sqrt(eval(n))*real(w((N*(k-1))+2:(N*(k-1))+1+N)'*evec(:,n));
    %             end
    %             s(m) = (1/sqrt(N_bar(m)))*t;
    %             sr(m) = sqrt(gamma)*norm(H{m,k}*w,2);%sum_square_abs(H{m,k}*w);
    %         end
    %     end
    % 
    %     minimize obj;
    %     subject to
    %     obj >= norm(w,2);%sum_square_abs(w);
    %     w(1) == 1;
    %     for m = 1:M
    %         s(m) >= sr(m);
    %     end
    %     cvx_end
    %     clear s sr t;
    % 
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
    %     if (cvx_optval > 0) && (ml == 10)
    %         [Mxm, MxI] = max(z);
    end

end

