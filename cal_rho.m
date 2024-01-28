function [rhoP1,rhoM1] = cal_rho(var,mean,pixel)
    zero = 10^(-100);
    wetCost = 10^8;
    beta = sqrt(var);
    shift = pixel-mean;
    divider_3 = normcdf((1.5+shift)./beta)-normcdf((-1.5+shift)./beta);
    pChangeP1_3 = (normcdf((1.5+shift)./beta)-normcdf((0.5+shift)./beta))./divider_3;
    pChangeM1_3 = (normcdf((-0.5+shift)./beta)-normcdf((-1.5+shift)./beta))./divider_3;
    pUnChange_3 = 1 - pChangeP1_3 - pChangeM1_3;
    rhoP1=log(pUnChange_3./(pChangeP1_3+zero));   %flip cost to probability
    rhoM1=log(pUnChange_3./(pChangeM1_3+zero)); 
    rhoM1(rhoM1>100)=wetCost;
    rhoP1(rhoP1>100)=wetCost;
end