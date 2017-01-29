% read the input digit in a specific location 
% function [i_r,j_r,k,diag,count,w,n_r,u_r,d_plus_in, d_minus_in] = D_in_control(d_plus,d_minus)
function [diag,count,k,u_r,n_r,q_plus_in, q_minus_in,q_o_plus,q_o_minus] = q_control_div(q_plus,q_minus)
unrolling = 8;
delta = 4;
persistent i; 
persistent j;
    if(isempty(i)&& isempty(j))
        i=0;  % left shift delta 0
        j=0;
    end
persistent diag_c;
persistent count_c;
persistent ite_k;
if isempty(diag_c) || isempty(count_c) || isempty(ite_k)
    diag_c = 0;
    count_c = 1;
    ite_k = 0;  % ite_k == 0, store the original x and y; ite_k > 0 store each iteration digits
end

persistent i_w; 
persistent j_w;
    if(isempty(i_w)&& isempty(j_w))
        % the reason i = -1 is to add the initial step in top function
        i_w=-1;  % left shift delta 0
        j_w=0;
    end
persistent d_plus_w;
persistent d_minus_w;
    if(isempty(d_plus_w)&& isempty(d_minus_w))
        d_plus_w = zeros(32,unrolling);  % left shift delta 0
        d_minus_w = zeros(32,unrolling);
        c = [1,0,1,0,0,1,1,0,1,1,1,1,1];
        d = [0,0,0,1,1,0,0,0,1,1,1,1,1];
        e=zeros(1,256);
        x_plus=[c,e];
        x_minus=[d,e];
        for o = 1:8*unrolling
            o_n = ceil(o/unrolling) - 1;
            if mod (o,unrolling)==0
                u_o = unrolling;
            else
                u_o = mod (o,unrolling);
            end
            d_plus_w(pairing(o_n,0),u_o) = x_plus(o);
            d_minus_w(pairing(o_n,0),u_o) = x_minus(o);
        end
    end
persistent diag_c_w;
persistent count_c_w;
persistent ite_k_w;
if isempty(diag_c_w) || isempty(count_c_w) || isempty(ite_k_w)
    diag_c_w = 0;
    count_c_w = 1;
    ite_k_w = 1;  % ite_k = 1 store 1st iteration digits
end

% digit output write
% in hardware count_c_w = 1,3,6,10,...to last digit, next clk diag_c + 1;
% in software,count_c_w = 1,2,4,7,...to first digit, current call diag_c + 1;  
if count_c_w == (diag_c_w+1)*diag_c_w/2 + 1
    diag_c_w = diag_c_w + 1;
end

if count_c_w == (diag_c_w-1)*diag_c_w/2 + 1
    ite_k_w = 1;   % generate 1st iteration digits 
else 
    % valid in 1st digit of delta group
    if i_w == 0 && j_w == 0
    ite_k_w = ite_k_w + 1;
    end
end

% digit assignment in delta group
if diag_c_w - ite_k_w == 0
    % import the first 4 digit in each iteration
    i_w = i_w + 1;
    if i_w == delta + 1
        count_c_w = count_c_w + 1;
        d_plus_w(pairing(0,ite_k_w),i_w-delta) = q_plus;
        d_minus_w(pairing(0,ite_k_w),i_w-delta) = q_minus;
        i_w = 0;
    end
    %d_o_plus = d_plus_w(pairing(0,ite_k_w),delta+1 : unrolling + delta);
    %d_o_minus = d_minus_w(pairing(0,ite_k_w),delta+1 : unrolling + delta);
    %d_o_plus = d_plus_w(pairing(0,ite_k_w),:);
    %d_o_minus = d_minus_w(pairing(0,ite_k_w),:);
    q_o_plus = d_plus_w;
    q_o_minus = d_minus_w;
    w_w = 0; n_w = 0;
else % diad_c - ite_k > 0
    % w_w is very important parameter to determine the digit location
    j_w = j_w+1;
        %w_w = (delta+1)+(diag_c_w-ite_k_w -1)*delta + j_w - delta;
        w_w = (diag_c_w-ite_k_w -1)*delta + j_w + 1;
        n_w = ceil(w_w/unrolling) - 1;
        if mod (w_w,unrolling)==0
            u_w = unrolling;
        else
            u_w = mod (w_w,unrolling);
        end
        %u_r = mod (w,unrolling);
        %d_minus_in = d_minus(ite_k,(delta+1)+(diag_c-ite_k)*delta + j);
        d_plus_w(pairing(n_w,ite_k_w),u_w) = q_plus;
        d_minus_w(pairing(n_w,ite_k_w),u_w) = q_minus;
    if j_w == delta
        count_c_w = count_c_w + 1;
        j_w = 0;
    end
    %d_o_plus = d_plus_w(pairing(n_w,ite_k_w),delta+1 : unrolling + delta);
    %d_o_minus = d_minus_w(pairing(n_w,ite_k_w),delta+1 : unrolling + delta);
    q_o_plus = d_plus_w;
    q_o_minus = d_minus_w;
end


% digit input read
% in hardware count_c = 1,3,6,10,...to last digit, next clk diag_c + 1;
% in software,count_c = 1,2,4,7,...to first digit, current call diag_c + 1;
if count_c == (diag_c+1)*diag_c/2 + 1
    diag_c = diag_c + 1;
end
    
if count_c == (diag_c-1)*diag_c/2 + 1    
    ite_k = 0;  % generate original input x,y 
else 
    % valid in 1st digit of delta group
    if i == 0 && j == 0
    ite_k = ite_k + 1;
    end
end

% digit assignment in delta group
if diag_c - ite_k == 1
    % import the first 4 digit in each iteration
    i = i + 1;   
    u_r=i;
        q_plus_in = d_plus_w(pairing(0,ite_k),i);
        q_minus_in = d_minus_w(pairing(0,ite_k),i);
    if i == delta + 1
        count_c = count_c + 1;
        i = 0;
    end
    w=0;n_r=0;
else % diad_c - ite_k > 1
    j = j+1;
        w = (delta+1)+(diag_c-ite_k - 2)*delta + j;
        n_r = ceil(w/unrolling) - 1;
        if mod (w,unrolling)==0
            u_r = unrolling;
        else
            u_r = mod (w,unrolling);
        end
        %u_r = mod (w,unrolling);
        %d_minus_in = d_minus(ite_k,(delta+1)+(diag_c-ite_k)*delta + j);
        q_plus_in = d_plus_w(pairing(n_r,ite_k),u_r);
        q_minus_in = d_minus_w(pairing(n_r,ite_k),u_r);
    %end
    if j == delta
        count_c = count_c + 1;
        j = 0;
    end

end
i_r=i;
j_r=j;
k=ite_k;
diag=diag_c;
count=count_c;
end



