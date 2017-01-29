% Digit-vector: x[P] = x_in(1,2,...,P); y[P] = y_in(1,2,...,P);
% x(j) = x(j)_plus - x(j)_minus; y(j) = y(j)_plus - y(j)_minus 
% if j > delta
% p = p(j-delta)
function [u,j_count,ena_q,CAd_plus,CAd_minus,CAq_plus_o,CAq_minus_o,CAq_plus_sel, CAq_minus_sel,CAd_plus_sel, CAd_minus_sel, shift_to_int_plus,shift_to_int_minus,v_frac_plus,v_frac_minus,compare_frac,w_frac_plus,w_frac_minus,shift_out_plus,shift_out_minus,count_one_plus,count_one_minus,q_plus,q_minus,v_int_plus,v_int_minus,w_int_plus,w_int_minus,q1_out,q0_out] = OD_FPT_ALG_v2(x_in_plus, x_in_minus, d_in_plus, d_in_minus, k)
% initialization step
% time_step: j
unrolling = 8;
delta = 4;
persistent j;
if isempty(j)
	j = 0;
end
persistent flag;
%persistent p_count;
% BRAM: w_frac
persistent w_plus_wr_frac;
persistent w_minus_wr_frac;
    if(isempty(w_plus_wr_frac)&& isempty(w_minus_wr_frac))
        w_plus_wr_frac=zeros(1024,unrolling);  
        w_minus_wr_frac=zeros(1024,unrolling);
    end
% BRAM: w_int
persistent w_plus_wr_int;
persistent w_minus_wr_int;
    if(isempty(w_plus_wr_int)&& isempty(w_minus_wr_int))
        w_plus_wr_int=zeros(1024,5);  
        w_minus_wr_int=zeros(1024,5);
    end
persistent q_plus_rev; persistent q_minus_rev;
persistent q_plus_comp; persistent q_minus_comp;
if isempty (q_plus_comp) || isempty(q_minus_comp)
    q_plus_comp = 0 ;%zeros(1024,1);
    q_minus_comp = 0; zeros(1024,1);
end
persistent CAw_frac_plus; persistent CAw_frac_minus;
persistent CAw_plus_int; persistent CAw_minus_int;
persistent q_out_plus; persistent q_out_minus;
% q[j_q] control signals
persistent j_q;
if isempty(j_q)
    j_q = 0;
end
persistent enable_q;
if isempty(enable_q)
    enable_q = 1;
end
persistent n_r_q;
persistent u_r_q;
persistent CAq_plus;
persistent CAq_minus;

if isempty(flag)
    %cin_one_plus = 0; cin_one_minus = 0; cin_two_plus = 0; cin_two_minus = 0; 
    CAw_frac_plus=zeros(1,unrolling); CAw_frac_minus=zeros(1,unrolling);
    q_out_plus = zeros(1,4*unrolling); q_out_minus = zeros(1,4*unrolling); 
    flag = 0;
end
                
for j = 1:64
        %j = j + 1;
        % x delay 3 clks in Aaron
		x_plus_rev = x_in_plus(j);
        x_minus_rev = x_in_minus(j);
        d_plus_rev = d_in_plus(j);
        d_minus_rev = d_in_minus(j);
        %enable = 1; % only valid when a new digit valid
        res_enable = 1;       
        n_r = ceil(j/unrolling) - 1;
        if mod(j,unrolling) == 0
            u_r = unrolling;
        else
            u_r = mod(j,unrolling);
        end
        for i = n_r:-1:0
            if i == n_r
                enable = 1;
            else
                enable = 0;
            end
            % initial delta step
            if (n_r==0)&&(u_r<=delta)
                % enable =1, write/read new digit, enable =0, read other chunk digits 
                % In OD, q[j] <=> x[j], y[j+1] <=> d[j+1]
                %[CAq_plus,CAq_minus,CAd_plus,CAd_minus] = CA_gen(q_plus_rev,q_minus_rev,d_plus_rev,d_minus_rev,i,i,u_r,k,enable);        
			    % CA[dj+1]
                [CAd_plus,CAd_minus] = CA_gen_d(d_plus_rev,d_minus_rev,i,i,u_r,k,enable);
                q_plus_rev = 0;
                q_minus_rev = 0;                
                CAq_plus_sel = zeros(1,8);
                CAq_minus_sel = zeros(1,8);   
			    [CAd_plus_sel, CAd_minus_sel]=SDVM_d(~q_plus_rev,~q_minus_rev,CAd_plus,CAd_minus);
                % Algorithm: v[j] =...; qj =...;
                % BRAM: w_frac read
                CAw_frac_plus = w_plus_wr_frac(pairing(i,k),:);
                CAw_frac_minus = w_minus_wr_frac(pairing(i,k),:);
                [shift_to_int_plus,shift_to_int_minus,v_frac_plus,v_frac_minus,compare_frac] = v_frac_div(CAq_plus_sel,CAq_minus_sel,CAw_frac_plus,CAw_frac_minus,x_plus_rev,x_minus_rev,  i,n_r);
                % BRAM: w_int read
                CAw_plus_int=w_plus_wr_int(pairing(i,k),:);
                CAw_minus_int=w_minus_wr_int(pairing(i,k),:);  
                [v_int_plus,v_int_minus] = v_int_div_no_q(compare_frac,CAw_plus_int,CAw_minus_int,shift_to_int_plus,shift_to_int_minus);
                % Algorithm: w[j+1]= v[j]; 
                % BRAM_frac write: 
                [w_frac_plus,w_frac_minus,shift_out_plus,shift_out_minus,count_one_plus,count_one_minus] = w_frac_div(CAd_plus_sel,CAd_minus_sel,v_frac_plus, v_frac_minus, i,i,n_r,k);
                %CAw_frac_plus = w_frac_plus;
                %CAw_frac_minus = w_frac_minus;
                w_plus_wr_frac(pairing(i,k),:)=w_frac_plus;
                w_minus_wr_frac(pairing(i,k),:)=w_frac_minus;
                % BRAM_int write
                [w_int_plus,w_int_minus] = w_int_div(CAd_plus_sel,CAd_minus_sel,v_int_plus, v_int_minus,i, shift_out_plus,shift_out_minus,count_one_plus,count_one_minus,k);
                w_plus_wr_int(pairing(i,k),:)=w_int_plus;
                w_minus_wr_int(pairing(i,k),:)=w_int_minus;     
                q_plus =0; q_minus=0;
                CAq_plus = zeros(1,unrolling); CAq_minus = zeros(1,unrolling);
            else
                % CA[dj+1]
                [CAd_plus,CAd_minus] = CA_gen_d(d_plus_rev,d_minus_rev,i,i,u_r,k,enable);
                
%                 % BRAM: q read
%                 q_plus_rev = q_plus_comp;
%                 q_minus_rev = q_minus_comp;
%                 % q[j_q] assignment
%                 j_q = unrolling * n_r + u_r - delta - 1;  
%                 n_q = ceil(j_q/unrolling) - 1;
%                 if mod(unrolling * n_r + u_r - delta - 1 , unrolling) == 0
%                     u_r_q = unrolling;
%                 else
%                     u_r_q = mod(unrolling * n_r + u_r - delta - 1 , unrolling);
%                 end 
%                 if n_q_i == n_q
%                     enable_q = 1;
%                 else
%                     enable_q = 0;
%                 end
%                 if j_q == 0 
%                     CAq_plus = zeros(1,unrolling); 
%                     CAq_minus = zeros(1,unrolling);
%                 else
%                     [CAq_plus,CAq_minus] = CA_gen_q(q_plus_rev,q_minus_rev,n_q_i,n_q_i,u_r_q,k,enable_q);
%                     if n_q_i >= 1
%                         n_q_i = n_q_i-1;
%                     else
%                         n_q_i = n_q;
%                     end
%                 end

                % Algorithm: CAx_sel = x[j]*y(j); CAy_sel = y[j+1]*x(j) 
                % in Aaron's verilog, there is a delta delay for CAq, may be ....  
			    [CAq_plus_sel, CAq_minus_sel]=SDVM_d(~d_plus_rev,~d_minus_rev,CAq_plus,CAq_minus);
			    %[CAd_plus_sel, CAd_minus_sel]=SDVM_d(~q_plus_rev,~q_minus_rev,CAd_plus,CAd_minus);
                      
                % Algorithm: v[j] =...; qj =...;
                % BRAM: w_frac read
                CAw_frac_plus = w_plus_wr_frac(pairing(i,k),:);
                CAw_frac_minus = w_minus_wr_frac(pairing(i,k),:);
                [shift_to_int_plus,shift_to_int_minus,v_frac_plus,v_frac_minus,compare_frac] = v_frac_div(CAq_plus_sel,CAq_minus_sel,CAw_frac_plus,CAw_frac_minus,x_plus_rev,x_minus_rev,  i,n_r);
                
                % BRAM: w_int read
                CAw_plus_int=w_plus_wr_int(pairing(i,k),:);
                CAw_minus_int=w_minus_wr_int(pairing(i,k),:);   
                [q_plus,q_minus,v_int_plus,v_int_minus] = v_int_div(compare_frac,CAw_plus_int,CAw_minus_int,shift_to_int_plus,shift_to_int_minus);
%                 % BRAM: q write
%                 q_plus_comp = q_plus;
%                 q_minus_comp = q_minus;      
                j_q = j - delta;
                j_count = j_q;
                ena_q = enable_q;
                if enable_q == 1
                    n_r_q = ceil(j_q/unrolling) - 1;
                    if mod(j_q,unrolling) == 0
                        u_r_q = unrolling;
                    else
                        u_r_q = mod(j_q,unrolling);
                    end
                    u = u_r_q;
                end                    
                [CAq_plus,CAq_minus] = CA_gen_q(q_plus,q_minus,n_r_q,n_r_q,u_r_q,k,enable_q);
                CAq_plus_o = CAq_plus; CAq_minus_o = CAq_minus;
                if n_r_q >= 1
                    n_r_q = n_r_q -1;
                    enable_q = 0;
                else
                    enable_q = 1;
                end          

                [CAd_plus_sel, CAd_minus_sel]=SDVM_d(~q_plus,~q_minus,CAd_plus,CAd_minus);
                % Algorithm: 2*w[j]=...; 
                % BRAM_frac write: 
                [w_frac_plus,w_frac_minus,shift_out_plus,shift_out_minus,count_one_plus,count_one_minus] = w_frac_div(CAd_plus_sel,CAd_minus_sel,v_frac_plus, v_frac_minus, i,i,n_r,k);
                %CAw_frac_plus = w_frac_plus;
                %CAw_frac_minus = w_frac_minus;
                w_plus_wr_frac(pairing(i,k),:)=w_frac_plus;
                w_minus_wr_frac(pairing(i,k),:)=w_frac_minus;
                % BRAM_int write
                [w_int_plus,w_int_minus] = w_int_div(CAd_plus_sel,CAd_minus_sel,v_int_plus, v_int_minus,i, shift_out_plus,shift_out_minus,count_one_plus,count_one_minus,k);
                w_plus_wr_int(pairing(i,k),:)=w_int_plus;
                w_minus_wr_int(pairing(i,k),:)=w_int_minus;
            end
                
                flag = flag + 1;
                count = j;
                q_out_plus(1,j)=q_plus;
                q_out_minus(1,j)=q_minus;
                q1_out=q_out_plus;
                q0_out=q_out_minus;

        end
        
end
%end
			

