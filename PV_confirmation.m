clear all
clc
load('pv_experiment.mat')
data=veriler;

% experimental SD model parameters for IRUN algorithm
% Iph=0.760787964314995;
% Isd=0.000000310690918;
% Rs=0.036546780767411;
% Rsh=52.889903680665235;
% n=1.477271627088029;

% experimental DD model parameters for IRUN algorithm
Iph=0.760802681986862;
Isd1=0.000001000000000;
Isd2=0.000000089070984;
Rs=0.037629585940514;
Rsh=55.999261921534078;
n1=1.823919149311149;
n2=1.380399049245109;

%%
q=1.60217646e-19;
k=1.38064852e-23;
T=33+273.15;
Vt=k*T/q; 
If_tahmin=zeros(size(data,2),1);
syms If_t
numberOfdiodes=2;  % enter for SD model: 1 and for DD model: 2
for i=1:size(data,1)
    if numberOfdiodes==1 
        eqn = If_t==Iph-Isd*(exp((data(i,1)+If_t*Rs)/(n*Vt))-1)-(data(i,1)+If_t*Rs)/Rsh;
    else if numberOfdiodes==2
        eqn = If_t==Iph-Isd1*(exp((data(i,1)+If_t*Rs)/(n1*Vt))-1)-Isd2*(exp((data(i,1)+If_t*Rs)/(n2*Vt))-1)-(data(i,1)+If_t*Rs)/Rsh;
        end
    end
     V = vpasolve(eqn,If_t);
     if size(V,1)==0
         V=100;
     end
     dbstop if error
     If_tahmin(i)=V;
end
e=If_tahmin-data(:,2);
iae=mean(abs(e))
rmse=sqrt(mean(e.^2))
[data(:,2) If_tahmin]