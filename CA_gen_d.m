%when call =j
%CAx1,CAx0 first read then write because q[j] will be used
%CAy1,CAy0 first write then read because d[j+1] will be used
% In OD, q[j] <=> x[j], y[j+1] <=> d[j+1]
function [CAd1,CAd0] = CA_gen_d(d1,d0,wr_addr,rd_addr,u_r,ite_input_r,enable)%,refresh)
unrolling = 8; 
persistent CA_d1;
persistent CA_d0;
    if(isempty(CA_d1)&& isempty(CA_d0))
        CA_d1=zeros(256*4,unrolling);  % 64*4
        CA_d0=zeros(256*4,unrolling);
    end
    if(enable == 1)
% refresh is not needed in our function coz every time we take a nre digit
% for computation.this case can be considered as "always refresh"
        %if(redresh==1 && rd_addr==0)    
% d, this step can be easily implemented with "Register", but in Matlab, it can be implemented as below     
            addr_r=pairing(rd_addr, ite_input_r);
            addr_w=pairing(wr_addr, ite_input_r);  
            CA_d1(addr_w,u_r) = d1;    %y[j+1];
            CA_d0(addr_w,u_r) = d0;
        %d,first write, then read    
            CAd1 = CA_d1(addr_r,1:unrolling);
            CAd0 = CA_d0(addr_r,1:unrolling);
    else
            addr_r=pairing(rd_addr,ite_input_r);
            CAd1 = CA_d1(addr_r,1:unrolling);
            CAd0 = CA_d0(addr_r,1:unrolling);
    end
end
        
            
