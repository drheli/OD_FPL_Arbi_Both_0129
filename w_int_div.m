function [w_plus_int,w_minus_int] = w_int_div(d_plus,d_minus,v_int_plus, v_int_minus,wr_addr, shift_out_plus,shift_out_minus,cout_one_plus,count_one_minus,k)
persistent shift_in_plus;
persistent shift_in_minus;
if isempty(shift_in_plus) || isempty(shift_in_minus)
    shift_in_plus = zeros(1,2);
    shift_in_minus = zeros(1,2);
end
persistent w_int_plus;
persistent w_int_minus;
    if(isempty(w_int_plus)&& isempty(w_int_minus))
        w_int_plus=zeros(256,5);  % 64*4
        w_int_minus=zeros(256,5);
    end
w_int_plus_rev = zeros(1,5);
w_int_minus_rev = zeros(1,5);
    
if wr_addr == 0
    shift_in_plus(2) = shift_out_plus; 
    shift_in_plus(1) = cout_one_plus;
    shift_in_minus(2) = shift_out_minus;
    shift_in_minus(1) = count_one_minus;
else
    shift_in_plus = zeros(1,2);
    shift_in_minus = zeros(1,2);
end

%if wr_enable ==1
    if wr_addr == 0
        temp_plus1 = dec2bin(bin2dec(num2str(v_int_plus)) + bin2dec(num2str(d_plus(1:3))))- '0'; 
    % w_int_plus_rev
    if length(temp_plus1)==6
        for i=1:5
            w_int_plus_rev(1,i) = temp_plus1(i+1);
        end
    end
    if length(temp_plus1)==5
        for i=1:5
            w_int_plus_rev(1,i) = temp_plus1(i);
        end
    end
    if length(temp_plus1)==4
        w_int_plus_rev(1,1) = 0;
        for i=1:4
            w_int_plus_rev(1,i+1) = temp_plus1(i);
        end
    end
    if length(temp_plus1)==3
        w_int_plus_rev(1,1:2) = zeros(1,2);
        for i=1:3
            w_int_plus_rev(1,i+2) = temp_plus1(i);
        end
    end
    if length(temp_plus1)==2
        w_int_plus_rev(1,1:3) = zeros(1,3);
        for i=1:2
            w_int_plus_rev(1,i+3) = temp_plus1(i);
        end
    end
    if length(temp_plus1)==1
        w_int_plus_rev(1,1:4) = zeros(1,4);
        w_int_plus_rev(1,5) = temp_plus1(1);
    end        
    
        temp_plus2 = dec2bin(bin2dec(num2str(w_int_plus_rev(2:5))) *2 + bin2dec(num2str(shift_in_plus))) - '0'; 
        %test = length(temp_plus2);
    % w_int_plus
    if length(temp_plus2)==6
        for i=1:5
            w_int_plus(pairing(wr_addr,k),i) = temp_plus2(i+1);
        end
    end
    if length(temp_plus2)==5
        for i=1:5
            w_int_plus(pairing(wr_addr,k),i) = temp_plus2(i);
        end
    end
    if length(temp_plus2)==4
        w_int_plus(pairing(wr_addr,k),1) = 0;
        for i=1:4
            w_int_plus(pairing(wr_addr,k),i+1) = temp_plus2(i);
        end
    end
    if length(temp_plus2)==3
        w_int_plus(pairing(wr_addr,k),1:2) = zeros(1,2);
        for i=1:3
            w_int_plus(pairing(wr_addr,k),i+2) = temp_plus2(i);
        end
    end
    if length(temp_plus2)==2
        w_int_plus(pairing(wr_addr,k),1:3) = zeros(1,3);
        for i=1:2
            w_int_plus(pairing(wr_addr,k),i+3) = temp_plus2(i);
        end
    end
    if length(temp_plus2)==1
        w_int_plus(pairing(wr_addr,k),1:4) = zeros(1,4);
        w_int_plus(pairing(wr_addr,k),5) = temp_plus2(1);
    end 
        
        temp_minus1 = dec2bin(bin2dec(num2str(v_int_minus)) + bin2dec(num2str(d_minus(1:3))))- '0'; 
    % w_int_plus_rev
    if length(temp_minus1)==6
        for i=1:5
            w_int_minus_rev(1,i) = temp_minus1(i+1);
        end
    end
    if length(temp_minus1)==5
        for i=1:5
            w_int_minus_rev(1,i) = temp_minus1(i);
        end
    end
    if length(temp_minus1)==4
        w_int_minus_rev(1,1) = 0;
        for i=1:4
            w_int_minus_rev(1,i+1) = temp_minus1(i);
        end
    end
    if length(temp_minus1)==3
        w_int_minus_rev(1,1:2) = zeros(1,2);
        for i=1:3
            w_int_minus_rev(1,i+2) = temp_minus1(i);
        end
    end
    if length(temp_minus1)==2
        w_int_minus_rev(1,1:3) = zeros(1,3);
        for i=1:2
            w_int_minus_rev(1,i+3) = temp_minus1(i);
        end
    end
    if length(temp_minus1)==1
        w_int_minus_rev(1,1:4) = zeros(1,4);
        w_int_minus_rev(1,5) = temp_minus1(1);
    end        
    
        temp_minus2 = dec2bin(bin2dec(num2str(w_int_minus_rev(2:5))) *2 + bin2dec(num2str(shift_in_minus)))- '0'; 
    % w_int_plus
    if length(temp_minus2)==6
        for i=1:5
            w_int_minus(pairing(wr_addr,k),i) = temp_minus2(i+1);
        end
    end
    if length(temp_minus2)==5
        for i=1:5
            w_int_minus(pairing(wr_addr,k),i) = temp_minus2(i);
        end
    end
    if length(temp_minus2)==4
        w_int_minus(pairing(wr_addr,k),1) = 0;
        for i=1:4
            w_int_minus(pairing(wr_addr,k),i+1) = temp_minus2(i);
        end
    end
    if length(temp_minus2)==3
        w_int_minus(pairing(wr_addr,k),1:2) = zeros(1,2);
        for i=1:3
            w_int_minus(pairing(wr_addr,k),i+2) = temp_minus2(i);
        end
    end
    if length(temp_minus2)==2
        w_int_minus(pairing(wr_addr,k),1:3) = zeros(1,3);
        for i=1:2
            w_int_minus(pairing(wr_addr,k),i+3) = temp_minus2(i);
        end
    end
    if length(temp_minus2)==1
        w_int_minus(pairing(wr_addr,k),1:4) = zeros(1,4);
        w_int_minus(pairing(wr_addr,k),5) = temp_minus2(1);
    end 
       
        w_plus_int = w_int_plus(pairing(wr_addr,k),:);
        w_minus_int = w_int_minus(pairing(wr_addr,k),:);
    else
        w_plus_int = zeros(1,5);
        w_minus_int = zeros(1,5);
    end
%end
        