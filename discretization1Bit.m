% Function of 1-bit quantization for RIS
function [Angle] = discretization1Bit(w)
Angle = angle(w);
Angle = Angle.*180/pi;
Angle(Angle>-180 & Angle<-90) = Angle(Angle>-180 & Angle<-90)+360;
Angle(Angle>-90 & Angle<=90) = 1;
Angle(Angle>90 & Angle<=270) = -1;
