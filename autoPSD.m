function [gkk] = autoPSD(w) 
global scal
   
wg=15.0;%natural frequency of the soil layer    %broad-band
ztg=0.6;%damping of the soil layer
wf=1.5;%frequency of the second filter
ztf=0.6;%damping of the second filter 

gkk=scal*(1+4.*ztg^2.*(w./wg).^2)./((1-(w./wg).^2).^2 + 4*ztg.^2*(w./wg).^2) .* ((w./wf).^4)./((1-(w./wf).^2).^2 + 4*ztf.^2*(w./wf).^2);
