function [velocity_now,bed,moving,mnQ]=transport400(velocity_now, bed, moving, i, n, p, xs, ys);
%this transport function is based on Bailard (1981)

%Copyright (C) <2007> <Edith Gallagher>  
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


%Sediment transport variables
D50=0.3; %mm
D50=D50./10; %mm -> cm
D90=0.5; %mm  %this is a total and complete guess
%goes pretty well with Truc Vert sand (if anything a little lower and
%narrower.  Obviously MRY and Sennen are more broadly distributed 7/20/07
D90=D90./10; %mm-> cm
g=980; %cm/sec**2
rhow=1.0;  %gm/cm**3
rhos=2.65;  %gm/cm**3
Cf=0.03;    %coeficient of friction %Cf=0.006%Cf=0.03
epsilonS=0.015; %unitless sus load efficiency factor
epsilonB=0.135; %unitless bedload efficiency factor
%tanphi=0.46;   %unitless angle of repose of 25 degree
tanphi=0.3;    %...for 17 degrees
nu=136.e-8;                           %viscosity
tanbeta=xs./10;     %bedslope x-only for now (m x n)
                    %this gives rise of xs over run of 10 cm
%
%  intermediate calculations 
%
Y=(rhos-rhow)*g;  %constant
%Dstr=((Y/nu^2)^0.33)*D50; %constant
%W=(nu/D50)*((107.33+1.049*Dstr^3)^0.5-10.36); %constant
W=3;
A=rhow*Cf*epsilonB/(Y*tanphi); %constant
%B=tanbeta./tanphi; %this is mxn
B=0;
C=rhow*Cf*epsilonS/(Y*W); %constant
%D=(epsilonS*tanbeta)/W; %this is mxn
D=0;

%further calculations
uraw=velocity_now;
%uraw(1:5,1:5)
%vraw=0;

uvec=uraw.*uraw;
%vvec=vraw.*vraw;
%vector=uvec+vvec; %
vector=uvec;

%not using bailard the way I used to, for time averaged flows
%now using Bailard's eqn 11 for instantaneous flow

term1=vector.*uraw;
%term1(1:5,1:5)
term2=abs(uraw).^3; %this is downslope so no sign
%term2(1:5,1:5)
term3=term2.*uraw;
%term3(1:5,1:5)
term4=abs(uraw).^5; %again downslope
%term4(1:5,1:5)

%%%%
%%%I think these slope terms (below) turn out to be not so important (2011)
%%%%
%Aug 16 2007
%I keep getting this wrong: the single q equation is needed for letting
%slope enhance or decrease transport (the signs work out)
%I might still come back and remove slope effects for small bedforms-later

% the slope term: I am going to say that if the bedforms are big (?) then
% the term will be important
%bed_rms=sqrt(mean(mean((bed-mean(mean(bed))).^2)));
%if bed_rms > 0;  %then use the slope term
   %but I need to get the sign right
   %vsign=sign(mean(mean(velocity_now(:,:))));
   %if vsign > 0 %then downstream slope is neg, the term is neg, 
                % but we want it to increase xprt so change sign
%       'pos velocity'
       q=A.*term1-A.*B.*term2+C.*term3-C.*D.*term4;
       %A.*term1(1:5,1:5)
       %A.*B(1:5,1:5).*term2(1:5,1:5)
       %C.*term3(1:5,1:5)
       %C.*D(1:5,1:5).*term4(1:5,1:5)
       
       %q(1:5,1:5)
       %pause
   %else %vsign < 0 then downstream slope is pos and will enhance xprt
%       'neg velocity'
    %   q=A.*term1+A.*B.*term2+C.*term3+C.*D.*term4;
   %end
%else %if bedforms are small then don't include the slope term at all
%    'small bedforms'
  %  q=A.*term1+C.*term3;
%end
   
     mnQ=mean(mean(q))
     mnU=mean(mean(velocity_now))
     
    
%     mn1=mean(mean(A.*term1));
%     mn2=mean(mean(A.*B.*term2));
%     mn3=mean(mean(C.*term3));
%     mn4=mean(mean(C.*D.*term4));
%     figure(15)
%     hold on
%     plot(i,mnQ,'k.')
%     hold on
%     plot(i,mnU,'b.')
%     plot(i,mn1,'r.')
%     plot(i,mn2,'g.')
%     plot(i,mn3,'c.')
%     plot(i,mn4,'m.')
    
    %NOW apply the calculated transport -> q=volume/(unit t*unit width)
    %  one of my blocks equals 10cm x 10 cm x 1cm high
    %  that is 100 cm^3 * 0.6 volume concentration of sand in stationary
    %  bed - IE 60 cm^3 of the block is sand
    %10/9/06 changing this to 10 cm x 10cm x 1mm = 10 cm^3
    %10 cm^3*0.6 =6
    numblocks=abs(round(q));
    %numblocks=numblocks./6;  %# of blocks to be move given by rounding to 
                                    %the nearest integer number
                                    
    nmmean=mean(mean(numblocks))
    %figure(15)
    %hold on
    %plot(i,nmmean.*10,'r*')
    
                                    
    %If there are some pickup and some drop locations in the grid we have
    %got to check them all
    
    %mass_conservation1=sum(sum(bed))+sum(sum(moving));
    %loss1=6250000-mass_conservation1;
    %meannumblocks=mean(mean(numblocks));
    
    for nnn=1:n
        for ppp=1:p
            if numblocks(nnn,ppp) < 1
                  %put everything down
                bed(nnn,ppp)=bed(nnn,ppp)+moving(nnn,ppp);
                moving(nnn,ppp)=moving(nnn,ppp)-moving(nnn,ppp);
            else
                moving(nnn,ppp)=moving(nnn,ppp)+numblocks(nnn,ppp);  %pickup the # of blocks 
                bed(nnn,ppp)=bed(nnn,ppp)-numblocks(nnn,ppp); %remove that sand from bed
            end
        end
    end
    
    %mass_conservation=sum(sum(bed))+sum(sum(moving));
    %loss=6250000-mass_conservation;

    %figure(10)
    %hold on
    %meanmov(i)=mean(mean(moving));
    %plot(i,meanmov(i),'cx')
    
   %except you can't erode a bed with no sand, at least not in this model!
   %so add sediment back into the bed matrix
    %ie you can't have negative numbers in bed matrix, just zero
    [nonnegbedindn,nonnegbedindp]=find(bed < 0);
    lennonnegbed=length(nonnegbedindn);
    if lennonnegbed > 0
        for a=1:lennonnegbed
            nnn=nonnegbedindn(a); %these two variable are just too darn long
            ppp=nonnegbedindp(a);
            moving(nnn,ppp)=moving(nnn,ppp)-abs(bed(nnn,ppp));
            bed(nnn,ppp)=bed(nnn,ppp)+abs(bed(nnn,ppp));
        end
    end
    %end of function