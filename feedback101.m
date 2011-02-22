function [velocity_now, velocity_next, bed]=feedback101(velocity_now,velocity_next,bed, i, xs, ys, n, p);
%this function uses the bed morphology to alter the velocity field for the
%next time step

%Copyright (C) <2008> <Edith Gallagher>  
%Developer can be contacted at 
%edith.gallagher@fandm.edu
%or
%Edith Gallagher
%Franklin and Marshall College
%PO Box 3003
%Lancaster, PA 17604-3003

%This program is free software; you can redistribute it and/or modify it under the terms of 
%the GNU General Public License as published by the Free Software Foundation; either 
%version 2 of the License, or (at your option) any later version.  
%This program is distributed in the hope that it will be useful, but WITHOUT ANY 
%WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%PARTICULAR PURPOSE. See the GNU General Public License for more details.  
%You should have received a copy of the GNU General Public License (gpl.txt) along 
%with this program; if not, go to http://www.gnu.org/licenses/gpl-2.0.txt or 
%write to the Free Software Foundation, Inc., 
%51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


%feedback100: changes velocity next according to bed morphology now 
%%%%FIRST%%%
%velocity on steep downstream slopes, like shadow zone
%behind a bedform experiencing flow separation. vel(t+1)=vel(t+1).\100
%%%%second%%%
%Increase velocity just
%downstream of this shadow zone velocity (t+1)=velocity(t+1)+ 0.1*velocity(t+1), 
%like a reattachement point with some turbulence.  (Doesn't do much)
%%%%%%THIRD%%%%%
%Upstream slopes experience some acceleration depending on height of
%bed above mean:  accel=(height of bed above mean) %%%/5 (or 1/2)
%vel(t+1)=vel(t+1)+accel

%FeedBack loop
clear FBindn FBindp
%let's let the bed effect the velocity
vel_next_orig=velocity_next(:,:);  %hold existing velocity array for next time point
                    %now adapt it
                    
%here we will let vel go to zero in the shadow of a bed form
%but we will also let the turbulence increase velocity slightly just
%downstream of the shadow zone
%5/8/2006  zero velocity is bad in model.  let's let velocity get small,
%say ./100

%xs and ys should be good from the last iteration through the smoothing loop
%get general vel direction - direction of oscillatory flow 
vsign=sign(mean(mean(velocity_now(:,:))));
if vsign > 0
    clear FBindn FBindn2 FBindp FBindp2
    %for positive velocity, find downstream slopes (they will be negative)
    [FBindn, FBindp]=find(xs < -2);
    %[FBindn, FBindp]=find(xs < -4);  %test for scale with 20x20cm blocks
    
    for fb=1:length(FBindn)
        %shadowshadowshadowshadow
        %decrease velocity in lee of bedform
        velocity_next(FBindn(fb), FBindp(fb))=vel_next_orig(FBindn(fb), FBindp(fb))./100;%put 0 on the upstream end of slope and ...
        shadowindn=FBindn(fb)+1; %slope at n is for (n+1)-n
        if shadowindn > n
            shadowindn=shadowindn-n; %periodic boundary
        end
        velocity_next(shadowindn,FBindp(fb))=vel_next_orig(shadowindn,FBindp(fb))./100; %put 0 on downstream end of slope too
        
        %wakewakewakewakewakereattachment (2009, this doesn't do much)
        %AND increase velocity slightly at reattachment point owing to turbulence
        turbindn=FBindn(fb)+2;
        if turbindn > n
            turbindn=turbindn-n; %periodic boundary
        end
        velocity_next(turbindn,FBindp(fb))=vel_next_orig(turbindn,FBindp(fb))+...
            vel_next_orig(turbindn,FBindp(fb)).*0.1;  %add 10% for turbulent wake
    end

    %now generate some acceleration on upstream positive slope
    [FBindn2, FBindp2]=find(xs > 0);  %positive slope
    %[FBindn2, FBindp2]=find(xs > 1);  %test for scale with 20x20cm blocks
    mnelev=mean(mean(bed));  %find mean elevation
    for fb=1:length(FBindn2)
        %for pos slope the uphill/downstream point
        downstreamaccelind=FBindn2(fb)+1;
        if downstreamaccelind > n
            downstreamaccelind=downstreamaccelind-n; %periodic boundary
        end
        accel2=(bed(downstreamaccelind,FBindp2(fb))-mnelev); %./2;%take a fraction of 
                                                        %height (above mean) and add it to the
                                                        %velocity  so 30cm
                                                        %amplitude bump
                                                        %gives velocity
                                                        %boost of 15cm/sec
        velocity_next(downstreamaccelind,FBindp2(fb))=vel_next_orig(downstreamaccelind,FBindp2(fb)) + accel2;
    end
    
else %for negative mean velocity
    clear FBindn FBindn2 FBindp FBindp2
    %for negative velocity, find downstream slopes (they will be positive)
    [FBindn, FBindp]=find(xs > 2);
    %[FBindn, FBindp]=find(xs > 4);  %test for scale with 20x20cm blocks

    %for fb=1:length(FBindn)
    for fb=length(FBindn):-1:1  %go backwards through indices so we are moving downstream & downslope
        %%shadowshadowshadowshadow
        %%decrease velocity in lee of bedform
        velocity_next(FBindn(fb),FBindp(fb))=vel_next_orig(FBindn(fb),FBindp(fb))./100; %put 0 on downstream end of slope and ...
        shadowindn=FBindn(fb)+1; %slope at n is for (n+1)-n
        if shadowindn > n
            shadowindn=shadowindn-n; %periodic boundary
        end
        velocity_next(shadowindn,FBindp(fb))=vel_next_orig(shadowindn,FBindp(fb))./100; %put 0 on the upstream end of slope too
   
        %%wakewakewakewakewakereattachment
        %%AND increase velocity slightly at reattachment point owing to turbulence
        turbindn=FBindn(fb)-1; %going downstream and downslope
        if turbindn < 1
            turbindn=n+turbindn; %periodic boundary
        end
        velocity_next(turbindn,FBindp(fb))=vel_next_orig(turbindn,FBindp(fb))+...
            vel_next_orig(turbindn,FBindp(fb)).*0.1;  %add 10% for turbulent wake
    end
    
    %now generate some acceleration on upstream negative slope
    %for neg slope the uphill/downstream point
    [FBindn2, FBindp2]=find(xs < 0);  %negative slope
    %[FBindn2, FBindp2]=find(xs < -1);  %test for scale with 20x20cm blocks
    mnelev=mean(mean(bed));  %find mean elevation
    for fb=1:length(FBindn2)
        accel=(bed(FBindn2(fb),FBindp2(fb))-mnelev); %./2;  %take a fraction of 
                                                        %height (above mean) and add it to the
                                                        %velocity  so 30cm
                                                        %amplitude bump
                                                        %gives velocity
                                                        %boost of 6cm/sec
        velocity_next(FBindn2(fb),FBindp2(fb))=vel_next_orig(FBindn2(fb),FBindp2(fb)) - accel;
    end
   
end