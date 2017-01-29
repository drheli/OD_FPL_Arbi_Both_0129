function [v_int_plus,v_int_minus] = v_int_div_no_q(compare_frac,w_int_plus,w_int_minus,shift_to_int_plus,shift_to_int_minus) %,wr_addr,rd_addr, ite_input_count)
    v_int_plus = zeros(1,5);
    v_int_minus = zeros(1,5);
    v_int_plus_rev = bin2dec(num2str(w_int_plus))+bin2dec(num2str(shift_to_int_plus));
    v_int_minus_rev = bin2dec(num2str(w_int_minus))+bin2dec(num2str(shift_to_int_minus));
    v_int_plus_reg = dec2bin(v_int_plus_rev)-'0';
    v_int_minus_reg = dec2bin(v_int_minus_rev)-'0';    
    %test = length(v_int_plus_reg);
% v_int_plus_reg  
    if length(v_int_plus_reg)==5
        for i=1:5
            v_int_plus(1,i) = v_int_plus_reg(i);
        end
    end
    if length(v_int_plus_reg)==4
        v_int_plus(1,1) = 0;
        for i=1:4
            v_int_plus(1,i+1) = v_int_plus_reg(i);
        end
    end
    if length(v_int_plus_reg)==3
        v_int_plus(1,1:2) = zeros(1,2);
        for i=1:3
            v_int_plus(1,i+2) = v_int_plus_reg(i);
        end
    end
    if length(v_int_plus_reg)==2
        v_int_plus(1,1:3) = zeros(1,3);
        for i=1:2
            v_int_plus(1,i+3) = v_int_plus_reg(i);
        end
    end
    if length(v_int_plus_reg)==1
        v_int_plus(1,1:4) = zeros(1,4);
        v_int_plus(1,5) = v_int_plus_reg(1);
    end
% v_int_minus_reg    
    if length(v_int_minus_reg)==5
        for i=1:5
            v_int_minus(1,i) = v_int_minus_reg(i);
        end
    end
    if length(v_int_minus_reg)==4
        v_int_minus(1,1) = 0;
        for i=1:4
            v_int_minus(1,i+1) = v_int_minus_reg(i);
        end
    end
    if length(v_int_minus_reg)==3
        v_int_minus(1,1:2) = zeros(1,2);
        for i=1:3
            v_int_minus(1,i+2) = v_int_minus_reg(i);
        end
    end
    if length(v_int_minus_reg)==2
        v_int_minus(1,1:3) = zeros(1,3);
        for i=1:2
            v_int_minus(1,i+3) = v_int_minus_reg(i);
        end
    end
    if length(v_int_minus_reg)==1
        v_int_minus(1,1:4) = zeros(1,4);
        v_int_minus(1,5) = v_int_minus_reg(1);
    end
    no_use=compare_frac;
            
end
    
    
%     v_sample = num2str(v_upper(1:3));
%     switch (v_sample)
%         case [0,1,1]
%             p_value1 = 1; p_value0 = 0;
%         case [0,1,0]
%             p_value1 = 1; p_value0 = 0;
%         case [0,0,1]
%             p_value1 = 1; p_value0 = 0;
%         case [0,0,0]
%             p_value1 = 0; p_value0 = 0;
%         case [1,1,1]
%             p_value1 = 0; p_value0 = 0;
%         case [1,1,0]
%             p_value1 = 0; p_value0 = 1;
%         case [1,0,1]
%             p_value1 = 0; p_value0 = 1;
%         case [1,0,0]
%             p_value1 = 0; p_value0 = 1;       
%     end