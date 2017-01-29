%when call =j
%CAx1,CAx0 first read then write because q[j] will be used
%CAy1,CAy0 first write then read because d[j+1] will be used
% In OD, q[j] <=> x[j], y[j+1] <=> d[j+1]
function [CAq1,CAq0,CAd1,CAd0] = CA_gen(q1,q0,d1,d0,wr_addr,rd_addr,u_r,ite_input_r,enable)%,refresh)
unrolling = 8;
% persistent ite_output_count;     %call function once, ite_output_count add one. FUNCTION as Counter
%     if(isempty(ite_output_count))
%         ite_output_count=1;
%     else
%         ite_output_count=ite_output_count+1;
%     end
% 
% persistent ite_count;
% persistent ite_input_count;  
%     if(isempty(ite_count))
%         ite_count=0;
%     end
%     if(isempty(ite_input_count))
%         ite_input_count=0;
%     end
%     if ite_output_count == (1 + (ite_count + 1)*ite_count / 2)  % 
% 		ite_count=ite_count+1;    % Diagonal Count, next one
%         ite_input_count = 0;      % Iteration Count  
% 	else
% 		ite_input_count = ite_input_count + 1;
%     end
% 	N_w=floor((ite_count-ite_input_count)/64);
% 	%u = 63 - (ite_count-N_depth*64) + ite_input_count;
%     u_w = (ite_count-N_w*64) - ite_input_count;        %begin from 1 to 64
% %CA_register    

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
% x, this step can be easily implemented with "Register", but in Matlab, it can be implemented as below  
% x,first read, then write.      
            addr_r=pairing(rd_addr,ite_input_r);
            CAq1 = CA_q1(addr_r,1:unrolling);
            CAq0 = CA_q0(addr_r,1:unrolling);
            addr_w=pairing(wr_addr, ite_input_r);
            CA_q1(addr_w,u_r) = q1;    %x[j+1]=>x[j]; initial x1=0,x2=xin1;
            CA_q0(addr_w,u_r) = q0;      
            CA_d1(addr_w,u_r) = d1;    %y[j+1];
            CA_d0(addr_w,u_r) = d0;
        %y,first write, then read    
            CAd1 = CA_d1(addr_r,1:unrolling);
            CAd0 = CA_d0(addr_r,1:unrolling);

    else
            addr_r=pairing(rd_addr,ite_input_r);
            CAq1 = CA_q1(addr_r,1:unrolling);
            CAq0 = CA_q0(addr_r,1:unrolling);
            CAd1 = CA_d1(addr_r,1:unrolling);
            CAd0 = CA_d0(addr_r,1:unrolling);
    end
        
%     %x[j+1]=>x[j]; initial x1=0,x2=xin1; only n =0;
%     if(enable==1)
%         if(u_r==1)
%             CAx1 = zeros(1,16);
%             CAx0 = zeros(1,16);
%         else
%             CAx1(1:(u_r-1)) = CAx1_rev(1:(u_r-1)); 
%             CAx1(u_r:16) = zeros(1,(16-u_r+1));
%             CAx0(1:(u_r-1)) = CAx0_rev(1:(u_r-1));
%             CAx0(u_r:16) = zeros(1,(16-u_r+1));
%         end
%     end
end
        
            
