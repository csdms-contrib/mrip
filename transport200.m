function [velocity_now,bed,moving]=transport200(velocity_now, bed, moving,i, A, T, n, p);
%transport 200 is based on Ribberink 1998 transport model

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
g=980; %cm/s^2
D50=0.2; %mm
D50=D50./10; %mm -> cm
D90=0.5; %mm  %this is a total and complete guess
D90=D90./10; %mm-> cm
rhos=2.65; %gm/cm^3
rhow=1.0;  %gm/cm^3
delta=(rhos-rhow)./rhow;

    %Ribberink 1998  RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
    
    %oscillatory flows  oooooooooooooooooooooooooooooooooooooo
    %fw=0.5;  %start with this value
    %ad=100;
    %difference=fw-0.3;
    %while difference>0.01  %iterate to get fw
    %    theta=(0.25.*rhow.*fw.*A.^2)./((rhos-rhow).*g.*D50);  %time-averaged magnitude of Shields parameter during the wave cycle
    %    ks=max([3.*D90,D50.*(1+6.*(theta-1))]);
    %    fwlast=fw;
    %if ks/ad < 0.63
    %    fw=exp((5.2.*(ks./ad).^0.194)-5.98);
    %elseif ks/ad >= 0.63
    %    fw=0.3;
    %end
    %difference=abs(fw-fwlast);  %check change in fw
    %end %while loop
    
    %fw=0.3;
    %fw=0.03;
    fw=0.01;
    
    theta_t=[0.5.*fw./(delta.*g.*D50)].*abs(velocity_now(:,:)).*velocity_now(:,:);
    theta_c=[0.5.*fw./(delta.*g.*D50)].*20.^2;  %20cm/sec is good for 10 s wave and 0.2 mm sand
    xprtsign(:,:)=theta_t(:,:)./abs(theta_t(:,:));
    theta_use=abs(theta_t(:,:))-theta_c;
    [temp1,temp2]=find(theta_use < 0); %because MATLAB makes complex numbers when this value is neg
    for tt=1:length(temp1)
        theta_use(temp1(tt),temp2(tt))=0;
    end
    phi(:,:)=11.*((theta_use).^1.65).*xprtsign(:,:);
    xprt_norm=sqrt(delta.*g.*(D50).^3);
    q(:,:)=phi(:,:).* xprt_norm; %q is in cm2/sec
    
    meanq=mean(mean(q))
    stdq=mean(std(q));
    meanu=mean(mean(velocity_now(:,:)))
    
    %figure(12)
    %hold on
    %plot(i,meanq,'k.')
    %hold on
    %plot(i,meanu,'.b')
   
    %steady flows  ssssssssssssssssssssssssssssssssssss
    %z0=ks/30;
    %z=20; %height of flow estimate above bed, another guesstimate
    %fc=2(0.4/ln(z/z0));
    %tb=0.5.*rhow.*fc.*velocity(:,:).^2;
    %theta=tb./((rhos-rhow).*g.*D50);
    %phi=m(theta=theta_c).^n;
    %q=phi.* xprt_norm;
    
    %combined flows  %cccccccccccccccccccccccccccccccccccc
    
    
    
    %NOW apply the calculated transport -> q=volume/(unit t*unit width)
    %  one of my blocks equals 10cm x 10 cm x 1cm high
    %  that is 100 cm^3 * 0.6 volume concentration of sand in stationary
    %  bed - IE 60 cm^3 of the block is sand
    %8/5/07 If the transport calculation gives cm^2/sec
    %numblocks=abs(round(q));
    %numblocks(:,:)=numblocks./60));  %# of blocks to be move given by rounding to 
                                    %the nearest integer number
    %8/6/07
    %now lets say that if the calculated transport is 1cm^3/sec the flow
    %picks up one block that is 1 cm^3 during this 1 sec time step.
    %so is it reasonable to say that is the flow acting at one (1 cm^2)
    %spot and so for that model unit 1cm^3/s picks up 1 block (10x10)?
    %Well that is what I am doing for now.
    numblocks=abs(round(q));
    
    %nmmean(i)=mean(mean(numblocks));
    %figure(12)
    %hold on
    %plot(i,nmmean(i).*10,'r*')
    
                                    
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