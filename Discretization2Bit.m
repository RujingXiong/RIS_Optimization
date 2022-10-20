function [Angle] = Discretization(w)
Angle = angle(w);
Angle = Angle.*180/pi;

Angle(Angle>360) = Angle(Angle>360)-360; %predeal
Angle(Angle<0) = Angle(Angle<0)+360;
Angle(Angle>360) = Angle(Angle>360)-360; %predeal
Angle(Angle<0) = Angle(Angle<0)+360;


Angle(Angle>315 | Angle<=45) = 1;%0
Angle(Angle>45 & Angle<=135) = 1i;
Angle(Angle>135 & Angle<=225) = -1;
Angle(Angle>225 & Angle<=315) = -1i;