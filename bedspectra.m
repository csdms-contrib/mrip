function [bed,spectra, plot_f]=bedspectra(bed, n, p, jump_frac, t,transport_mod,i);
%calculate the spectrum of the bedforms, plot it, and save it to a file

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


%calculate length scales with fft for every time step
%watch spectra evolve in time
%note this is Mark Orzech's code because Matlab's was wrong!
freq=0.1; %1 measurement every 10 cm
nfft=256;
ntot=n;
%nsegs=m;
nyq=(nfft/2)+1;
win=n;
hann=hanning(win);
dt=1/freq;
df=1/(nfft*dt);
plot_f=1:(nfft/2)+1;
plot_f=(plot_f*df)-(df); %this puts first spec est in 0 Hz freq bin
offset=0;

for s=1:p

ztemp1=bed(:,s)';
ztemp=detrend(ztemp1);

X=ztemp'.*hann;
ffttemp=fft(X,nfft);
Gx1(:,s)=((abs(ffttemp/nfft)).^2).*2.7/df; 
%2.7=correction factor for hanning window
%Gx(nyq+1:nfft,s)=[];
Gx(2:nyq,s)=2*Gx1(2:nyq,s);

end  %s loop on columns

ave_Gx=mean(Gx');

%figure(3)
%hold on
%semilogy(plot_f, ave_Gx,'b')
%axis([0 0.01 10000 3000000])


%write resulting data to a file
spectra=[ave_Gx'];
%pause

fid=fopen(['mrip_tprt' num2str(transport_mod) '_t' num2str(t) '_A' num2str(A) 'S' num2str(S) 'T' num2str(T) '_jf' num2str(jump_frac) '_spectra'],'a');
fwrite(fid,spectra,'double');
fclose(fid);

if i == t-1 %if we are in the last iteration then save x axis wavenumbers to the spectra file
    spec=[spectra plot_f'];
    fid=fopen(['mrip_tprt' num2str(transport_mod) '_t' num2str(t) '_A' num2str(A) 'S' num2str(S) 'T' num2str(T) '_jf' num2str(jump_frac) '_spectra'],'a');
    fwrite(fid,spec,'double');
    fclose(fid);
end

