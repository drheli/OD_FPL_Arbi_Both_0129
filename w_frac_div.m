function [w_frac_plus,w_frac_minus,shift_out_plus,shift_out_minus,count_one_plus,count_one_minus] = w_frac_div(d_plus,d_minus,v_frac_plus, v_frac_minus, wr_addr,rd_addr,n_r,k)
persistent cin_one_plus;
persistent cin_one_minus;
persistent shift_plus;
persistent shift_minus;
if isempty(shift_plus) || isempty(shift_minus)
    shift_plus = 0;
    shift_minus = 0;
end
unrolling = 8;
persistent residue1;
persistent residue0;
    if(isempty(residue1)&& isempty(residue0))
        residue1=zeros(256,unrolling);  % 64*4
        residue0=zeros(256,unrolling);
    end

if rd_addr == n_r  % n > 0
    cin_one_plus = 0;  
    cin_one_minus = 0;
end

if rd_addr == 0
    d_plus_append = zeros(1,3);
    d_minus_append = zeros(1,3);
else 
    d_plus_append = d_plus(1:3);
    d_minus_append = d_minus(1:3);
end
d_plus_comp = [d_plus(4:8),d_plus_append];
d_minus_comp = [d_minus(4:8),d_minus_append];
[w_frac_plus_store,w_frac_minus_store,count_one_plus,count_one_minus] = fourbitadder(d_plus_comp, d_minus_comp, v_frac_plus, v_frac_minus, cin_one_plus, cin_one_minus);
    
if rd_addr >=1  % n > 0
    cin_one_plus = count_one_plus;  
    cin_one_minus = count_one_minus;
else
    cin_one_plus = 0;
    cin_one_minus = 0;
end

% wr_enable is very complicated in Aaron
%if wr_enable == 1
    if wr_addr == n_r  % different shift_ena =1(i=n_r) in OM
        shift_in_plus = 0;
        shift_in_minus = 0;
    else
        shift_in_plus = shift_plus;
        shift_in_minus = shift_minus;    
    end
    residue1(pairing(wr_addr,k),1:unrolling-1) = w_frac_plus_store(2:unrolling);
    residue1(pairing(wr_addr,k),unrolling) = shift_in_plus;    
    residue1(pairing(wr_addr,k),1:unrolling-1) = w_frac_minus_store(2:unrolling);
    residue1(pairing(wr_addr,k),unrolling) = shift_in_minus; 
    w_frac_plus = residue1(pairing(rd_addr,k),:);
    w_frac_minus = residue1(pairing(rd_addr,k),:);
%end
shift_plus = w_frac_plus_store(1);
shift_minus = w_frac_minus_store(1);
shift_out_plus = w_frac_plus_store(1);
shift_out_minus = w_frac_minus_store(1);
end