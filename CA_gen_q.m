%when call =j
%CAx1,CAx0 first read then write because q[j] will be used
%CAy1,CAy0 first write then read because d[j+1] will be used
% In OD, q[j] <=> x[j], y[j+1] <=> d[j+1]
function [CAq1,CAq0] = CA_gen_q(q1,q0,wr_addr,rd_addr,u_r_q,ite_input_r,enable_q)%,refresh)
unrolling = 8; 
persistent CA_q1;
persistent CA_q0;
% persistent CA_x1_sel;
% persistent CA_x0_sel;
    if(isempty(CA_q1)&&isempty(CA_q0))
        CA_q1=zeros(256*4,unrolling);  % 64*4
        CA_q0=zeros(256*4,unrolling);
%         CA_x1_sel=zeros(256,17);  
%         CA_x0_sel=zeros(256,17);
    end
    if(enable_q == 1)
% refresh is not needed in our function coz every time we take a nre digit
% for computation.this case can be considered as "always refresh"
        %if(redresh==1 && rd_addr==0)    
% x, this step can be easily implemented with "Register", but in Matlab, it can be implemented as below  
% q,first read, then write.      
            addr_r=pairing(rd_addr,ite_input_r);
            CAq1 = CA_q1(addr_r,1:unrolling);
            CAq0 = CA_q0(addr_r,1:unrolling);
            addr_w=pairing(wr_addr, ite_input_r);
            CA_q1(addr_w,u_r_q) = q1;    %x[j+1]=>x[j]; initial x1=0,x2=xin1;
            CA_q0(addr_w,u_r_q) = q0;      

    else
            addr_r=pairing(rd_addr,ite_input_r);
            CAq1 = CA_q1(addr_r,1:unrolling);
            CAq0 = CA_q0(addr_r,1:unrolling);
    end
end
        
            
