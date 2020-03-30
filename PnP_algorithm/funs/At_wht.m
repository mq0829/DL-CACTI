function  At = At_wht(z, P,Mea_num,row, col)
   
 len = row*col; 
 z_all = [z; zeros(len-Mea_num,1)];
 At = myfwht(z_all);
 At= At(P)./(row*col);
% A = ifwht(z_all);
 %A = temp(1:Mea_num);
 
end