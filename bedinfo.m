function [bed,moving, data]=bedinfo(bed, moving, A, S, T, var_amp, i, jump_frac,t,transport_mod);
%calculate mass conservation, bed rms and the like

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

%pause
%mass_conservation=sum(sum(bed))+sum(sum(moving))
%loss=250000-mass_conservation

%bed RMS
%bed_demean=bed-mean(mean(bed));
%bed2=bed_demean.^2;
%bed2_mean=mean(mean(bed2));
bed_rms=sqrt(mean(mean((bed-mean(mean(bed))).^2)));
data=[var_amp A S T jump_frac bed_rms i];


%fid=fopen(['testing_info'],'a');
fid=fopen(['mrip_tprt' num2str(transport_mod) '_t' num2str(t) '_A' num2str(A) 'S' num2str(S) 'T' num2str(T) 'jf' num2str(jump_frac) '_info'],'a');
fwrite(fid,data,'double');
fclose(fid);


