function [bed]=plotbed(bed,velocity, moving, xplot, yplot, i);
%plots each new bed with time iteration 

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
    
bed_plot=bed;   %*0.1; %1 mm per bed unit
bed_plot=bed_plot-mean(mean(bed_plot)); %fluc about mean bed
figure(10)
clf
%axes('position',[0.2 0.3 0.37 0.3]) %trying to make it small for RCEM
set(gca,'FontSize',12)
surf(yplot,xplot,bed_plot','LineStyle','none') %.5 cm per bed unit %fluc about mean bed
axis([0.1 max(yplot) 0.1 max(xplot)])
caxis([-40 40])
set(gca,'ZLim',[-40 40]);
%view(-37.5,50)
view(0,90)
%ylabel('Cross-flow axis (m)','FontSize',16)
%xlabel('Alongflow axis (m)','FontSize',16)
ylabel('Cross-flow axis (m)','FontSize',12)
xlabel('Alongflow axis (m)','FontSize',12)
colorbar('NorthOutside','FontSize',12)
text(10,31,'Bed Elevation (cm)','FontSize',12) %for RCEM
%colorbar('FontSize',12)
%text(27, 5, 'Bed Elevation (cm)', 'rotation',90,'FontSize',16) %for a regular plot
%title('turbulence amplitude = 0 cm/sec','FontSize',18)
box on

%F(i)=getframe;
clear bed_plot




