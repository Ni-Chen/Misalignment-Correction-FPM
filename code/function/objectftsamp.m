function objftsamp1= objectftsamp(objectRecoverFT,ss,m1)
kxc=ss(1);kyc=ss(2); n1=m1;

 kyl = round(kyc-(m1-1)/2);kyh = round(kyc+(m1-1)/2);
 kxl = round(kxc-(n1-1)/2);kxh = round(kxc+(n1-1)/2);

 objftsamp1=objectRecoverFT(kyl:kyh,kxl:kxh);
        
end







        
        
   

















