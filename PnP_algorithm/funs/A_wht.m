function  A = A_wht(z, Q,Mea_num)
   
temp = myfwht(z(Q));
% temp =fwht(z(:));
 A = temp(1:Mea_num);
 
end