%% GCI Calculator

format long

r = 2;
% f1 =0.373372713685821
% f2 =0.373375677478581
% f3 = 0.373411752534933
FSEC = 1.25;
f3= 9.86798754568237
f2= 9.91136520322797
f1= 9.93797586818314


p = log(abs(f3-f2)/abs(f2-f1))/log(r)
f_h0 = f1+((f1-f2)/(r^p -1))
GCI_12 = FSEC *(abs((f1-f2)/f1)/(r^p -1))
GCI_23 = FSEC *(abs((f2-f3)/f2)/(r^p -1))

AsymptoticConvergence = (GCI_12/GCI_23) * r^p

% %Final answer is fh=0 +/- GCI12 in percent
% %time and k variables previously stores
% plot(TimeValues, k_values);
% hold on
% plot(TimeValues_96,k_values_96)
% plot(TimeValues_48,k_values_48)
% 
% set(gca,'fontsize',16);
% xlabel('time'); ylabel('k(t)');
% title(' k(t) vs t');
% 
% legend({'192x128', '96x64', '48x32'}, Location='best')