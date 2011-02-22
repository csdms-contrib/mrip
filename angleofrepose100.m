function [bed, xs, ys]=angleofrepose100(bed,moving, n, p, yplot, xplot, i);

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

%testing for angle of repose
    pile_high=3; %cm this is for a 10cm x 10cm block so the slope is 3/10, 
                 %which gives an angle of 17 degrees
    %pile_high=6;    %testing for a scale change : 20x20cm block is 6/20 so the same
    
    %first put first row and column or bed in last position ie (n+1)x(p+1)
    %so when diff works we will get an nxp array and it will be periodic in
    %slope
    
    %NOTE got my p's, n's X's and Y's mixed up.  I am only changing the way
    %that I take "diff"  SO now
    %   +X
    %   +n
    %    |     this is confusing!
    %    |
    %    |
    %    \ /--------> +p, +Y
    
    %b_old_temp=bed(:,25)'
    
    
    bed_temp(1:n,1:p)=bed;
    bed_temp(n+1,1:p)=bed(1,:);
    bed_temp(1:n,p+1)=bed(:,1);
    bed_temp(n+1,p+1)=bed(n,p);
    %calculate slopes
    Xslope=diff(bed_temp);
    xs=Xslope(:,1:p);
    Yslope=diff(bed_temp')'; %giving the slope in the p or Y direction
    ys=Yslope(1:n,:);
    %Now check to see if we need to smooth
    count=0; %to keep track of how often we go through smoothing loop
             %I think it should go a couple times then be smooth
    all_slopes=[xs;ys];  %lump the slopes together
    clear smooth_ind
    smooth_ind=find(abs(all_slopes) > pile_high); %check for big slopes
    gothroughloop=isempty(smooth_ind);
    %gothroughloop=1;  %skip smoothing loop
    
while gothroughloop == 0  %if there are big slopes go through the smoothing loop
    i;
    %message=['in the smoothing loop']
    %mass_conservation=sum(sum(bed))+sum(sum(moving))
    %loss=225000000-mass_conservation
    %pause
    
    %use most recent bed
    %make it bigger by one unit for the diff routine
    bed_temp(1:n,1:p)=bed;
    bed_temp(n+1,1:p)=bed(1,:);
    bed_temp(1:n,p+1)=bed(:,1);
    bed_temp(n+1,p+1)=bed(n,p);
    %calculate slopes
    Xslope=diff(bed_temp);
    %now make it smaller again only n x p
    xs=Xslope(:,1:p);
    [xsindn,xsindp]=find(abs(xs) > pile_high);
    
    %go randomly through the indices instead of from 1:size(indn)
    nn=size(xsindn);
    pp=randperm(nn(1));  %rearrange the indices randomly
    for jj=1:size(pp')
        j=pp(jj);
        %move a block from top of steep hill to bottom
        xdir=sign(xs(xsindn(j),xsindp(j)));  %gives plus or minus one
        bed(xsindn(j),xsindp(j))=bed(xsindn(j),xsindp(j))+xdir;
        if xsindn(j)+1 > n
            useind=1; %periodic boundary
        else
            useind=xsindn(j)+1;
        end
        bed(useind,xsindp(j))=bed(useind,xsindp(j))-xdir;
    end;
    
    %comafterx=sum(sum(bed))+sum(sum(moving))
    %pause
    
    %use most recent bed
    bed_temp(1:n,1:p)=bed;
    bed_temp(n+1,1:p)=bed(1,:);
    bed_temp(1:n,p+1)=bed(:,1);
    bed_temp(n+1,p+1)=bed(n,p);
    %calculate Y slopes
    Yslope=diff(bed_temp')'; %giving the slope in the p or Y direction
    ys=Yslope(1:n,:);
    [ysindn,ysindp]=find(abs(ys) > pile_high);
    
    %go randomly through the indices instead of from 1:size(indn)
    nn=size(ysindn);
    pp=randperm(nn(1));  %rearrange the indices randomly
    for jj=1:size(pp')
        j=pp(jj);
        %move a block from top of steep hill to bottom
        ydir=sign(ys(ysindn(j),ysindp(j)));
        bed(ysindn(j),ysindp(j))=bed(ysindn(j),ysindp(j))+ydir;
        if ysindp(j)+1 > p
            useind=1; %periodic boundary
        else
            useind=ysindp(j)+1;
        end
        bed(ysindn(j),useind)=bed(ysindn(j),useind)-ydir;
    end
    
    %bed_after_both=bed(:,25:28)
    
    %comaftery=sum(sum(bed))+sum(sum(moving))
    %pause
                             
    count=count+1;
    clear xsindn xsindp ysindn ysindp xdir ydir
    clear xdir ydir xs ys Xslope Yslope
    clear all_slopes smooth_ind
    
    %use most recent bed, do we have to go through the loop again?
    bed_temp(1:n,1:p)=bed;
    bed_temp(n+1,1:p)=bed(1,:);
    bed_temp(1:n,p+1)=bed(:,1);
    bed_temp(n+1,p+1)=bed(n,p);
    %calculate slopes
    Xslope=diff(bed_temp);
    xs=Xslope(:,1:p);
    Yslope=diff(bed_temp')'; %giving the slope in the p or Y direction
    ys=Yslope(1:n,:);
    %Now check to see if we need to smooth
    all_slopes=[xs;ys];  %lump the slopes together
    smooth_ind=find(abs(all_slopes) > pile_high); %check for big slopes
    %size(smooth_ind)
    gothroughloop=isempty(smooth_ind);
    
    
end %while gothroughloop loop on whether there are big slopes or not
%count

