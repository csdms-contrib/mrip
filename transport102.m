function [velocity,bed,moving]=transport102(velocity, bed, moving,i);
%tansport 100 is the simplest transport routine - 5 velocity catagories for
%picking up or dropping sand

%Copyright (C) <2010> <Edith Gallagher>  
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


    %establish velocity ranges in which sediment is either picked up or
    %dropped.  Higher velocities pick up more sediment than medium
    %velocities, low velocities drop sediment
    
    %July 2010 these thresholds were change in response to reviewers (they
    %asked what they were based on, I did more back-o-the-envelope
    %calculations and changed them a little
    %This is all for W=1.5cm/sec which is for about 200 micron sand
    [indn_alldown,indp_alldown]=find(abs(velocity(:,:))<=50);  %high settling rates
    len_alldown=length(indn_alldown);
    [indn_onedown,indp_onedown]=find(abs(velocity(:,:))>50 & abs(velocity(:,:))<=70); %settling
    len_onedown=length(indn_onedown);
    [indn_oneup,indp_oneup]=find(abs(velocity(:,:))>70 & abs(velocity(:,:))<=85); %erosion 1
    len_oneup=length(indn_oneup);
    [indn_twoup,indp_twoup]=find(abs(velocity(:,:))>85 & abs(velocity(:,:))<=100); %erosion 2
    len_twoup=length(indn_twoup);
    [indn_3up, indp_3up]=find(abs(velocity(:,:))>100);  %erosion 3
    len_3up=length(indn_3up);
    
    %when veolcity is low, sediment settles, at the present location
   if len_alldown > 0
   for a=1:len_alldown
       bed(indn_alldown(a),indp_alldown(a))=bed(indn_alldown(a),indp_alldown(a))+moving(indn_alldown(a),indp_alldown(a));
       moving(indn_alldown(a),indp_alldown(a))=moving(indn_alldown(a),indp_alldown(a))-moving(indn_alldown(a),indp_alldown(a));  %put it all down
   end
   end
   if len_onedown > 0
   for a=1:len_onedown
    %deposit one block
    %moving(indn_onedown(a),indp_onedown(a))=moving(indn_onedown(a),indp_onedown(a))-1;
    %bed(indn_onedown(a),indp_onedown(a))=bed(indn_onedown(a),indp_onedown(a))+1;
    %deposit half of what is in water colum
    bed(indn_onedown(a),indp_onedown(a))=bed(indn_onedown(a),indp_onedown(a))+fix(moving(indn_onedown(a),indp_onedown(a))./2);
    moving(indn_onedown(a),indp_onedown(a))=moving(indn_onedown(a),indp_onedown(a))-fix(moving(indn_onedown(a),indp_onedown(a))./2);
    
   end
   end
    
   %when velocity is high, the bed erodes
    %so sediment is taken from bed matrix to the moving water matrix
   if len_oneup > 0
   for a=1:len_oneup
    moving(indn_oneup(a),indp_oneup(a))=moving(indn_oneup(a),indp_oneup(a))+1;
    bed(indn_oneup(a),indp_oneup(a))=bed(indn_oneup(a),indp_oneup(a))-1;
   end
   end
   if len_twoup > 0
   for a=1:len_twoup
    moving(indn_twoup(a),indp_twoup(a))=moving(indn_twoup(a),indp_twoup(a))+2;
    bed(indn_twoup(a),indp_twoup(a))=bed(indn_twoup(a),indp_twoup(a))-2;
   end
   end
   if len_3up > 0
   for a=1:len_3up
    moving(indn_3up(a),indp_3up(a))=moving(indn_3up(a),indp_3up(a))+3;
    bed(indn_3up(a),indp_3up(a))=bed(indn_3up(a),indp_3up(a))-3;
   end
   end
   
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