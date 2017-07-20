function [trueSlope1,trueSlope2,confSlope1,confSlope2,pvalSlope,stdSlope1,stdSlope2]=bootStrapDemingConfPNAS(vals1_1,vals1_2,vals2_1,vals2_2,iters)
%%%% calculates the correlations' slope using deming(total least squares) regression, and the bootstrap 
[btsr1,bootsam] = bootstrp(iters,@deming,vals1_1',vals1_2');
[btsr2,bootsam] = bootstrp(iters,@deming,vals2_1',vals2_2');
slope1=btsr1(:,2);
slope2=btsr2(:,2);
[hhh pvalSlope]=ttest(slope1,slope2);
stdSlope1=std(slope1);stdSlope2=std(slope2);
  [b_dem sigma x_est y_est stats]  = deming(vals1_1',vals1_2');
  trueB1=b_dem(1);
   trueSlope1=b_dem(2);
  [b_dem sigma x_est y_est stats]  = deming(vals2_1',vals2_2');
  trueB2=b_dem(1);
  trueSlope2=b_dem(2);
  
  confSlope1=prctile(slope1,[05 95])-trueSlope1;
  confSlope2=prctile(slope2,[05 95])-trueSlope2;
