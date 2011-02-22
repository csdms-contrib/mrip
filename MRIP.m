%mrip_model101  A self-organization model to generate nearshore megaripples  

%Copyright (C) <2011> <Edith Gallagher>  
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

%MRIP -  this model consists of three layers. The first layer represents
%the stationary bed. From this layer, sand blocks are picked up or put down
%to change the elevation at each location (m,n) in the matrix. The second
%layer represents the sand that has been lifted from the bed and is in
%motion. The third layer represents the flow field that drive the motion of
%the sand. The flow field  is prescribed by the user and is the same at all
%locations for a given time point, except for a random spatial fluctuation.

%INITIAL PARAMETERS

%CREATE ARRAYS
%m=100; %z or depth of sand bed
n=256; %x, cross-shore flow (the oscillatory wave motions are back and forth along this axis)
p=256; %y, alongshore flow (perpendicular to the wave flow)
t=800;  %time (number of steps in the simulation, equal to seconds at this time)
%t=86400; %one whole day

%create stationary bed array
bed=repmat(10000,[n p]);
%load A75S20T800SR_CSDMS %or load a mat file with an existing bed morphology
%bedvar=rand(n,p); %6/9/08 add a random fluctuation to the bed instead of the flow (reviewer suggestion)
%bedvar=round(bedvar);
%bed=bed+bedvar.*3;
bed_temp=ones([n+1,p+1]); %for calculating slopes
%get starting slope ssssssssssssssss
bed_temp(1:n,1:p)=bed;
bed_temp(n+1,1:p)=bed(1,:);
bed_temp(1:n,p+1)=bed(:,1);
bed_temp(n+1,p+1)=bed(n,p);
%calculate slopes
Xslope=diff(bed_temp);
xs=Xslope(:,1:p);
Yslope=diff(bed_temp')'; %giving the slope in the p or Y direction
ys=Yslope(1:n,:);

%ssssssssssssssssssssssssssssss
moving=zeros([n,p]);  %this will hold all the particles that have been lifted from the stationary bed (that are moving)
jump=zeros([n,p]);   %this will accumulate a distance for the above particles to move

%arrays for plotting
xplot=[1:n]*10; %10 cm per bed unit
xplot=xplot/100; %m
yplot=[1:n]*10; %10 cm per bed unit
yplot=yplot/100; %m

%flow parameters
A=75; %cm/sec oscillatory velocity amplitude
T=10; %sec wave period
w=2*pi/T; %wave frequency
S=20; %steady current magnitude (aligned with waves)

var_amp=[30]; %DOUBLE the trough to crest magnitude of random turbulent vel 
%fluctuation added to the velocity signal (see below, it is multiplied by -0.5 - 0.5)

%Jump frac is the distance the sand moves forward. It is a fraction of the
%velocity. %10 is minimum, see movesandinsus for argument
jump_frac=20;
    
vel_strt=ones([n p]);
    
%we need a random spatial variation in the flow to start the bedform building process
vel_var=rand(n,p);
vel_var=vel_var-0.5;  %make it go from -.5 to +.5
vel_var_next=rand(n,p);
vel_var_next=vel_var_next-0.5;  %make it go from -.5 to +.5

%velocity_now holds the velocity at t
velocity_now(:,:)=vel_strt(:,:).*A.*sin(w.*1); %+vel_strt(:,:).*Aig.*sin(wig.*1); %IG
velocity_now(:,:)=velocity_now(:,:)+(vel_var(:,:).*var_amp)+S; %co-directional
%velocity_next holds the velocity at t=t+1 (so we can alter it with feedback loop)
velocity_next(:,:)=vel_strt(:,:).*A.*sin(w.*(2)); %+vel_strt(:,:).*Aig.*sin(wig.*2); %IG
velocity_next(:,:)=velocity_next(:,:)+(vel_var_next(:,:).*var_amp)+S; %co-directional

rms=sqrt(mean(velocity_now.^2)); %just to keep track of stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The time step loop: Run the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totcnt=0;
cnt=0;
for i=1:t-1
    
    if i > 1
        vel_var_next=rand(n,p); %a new random vel fluctuation
        vel_var_next=vel_var_next-0.5;  %make it go from -.5 to +.5
        velocity_now=velocity_next;  %use modified velocity, after the feedback routine.
        velocity_next(:,:)=vel_strt(:,:).*A.*sin(w.*(i+1));
        velocity_next(:,:)=velocity_next(:,:)+(vel_var_next(:,:).*var_amp)+S; %co-directional
    end
    
    %figure(1)
    %hold on
    %plot(i,mean(mean(velocity_now(:,:))),'.')
    
    [i A S var_amp jump_frac] %print some info to the screen
    
    %%%%%%%%%%%%%
    %Call the transporting subroutine
    %this routine decides whether and how much sand to lift into suspension
    %%%%%%%%%%%%%
    %Crude transport model
    %New simple rules July 2010
    %[velocity_now,bed,moving]=transport102(velocity_now,bed,moving,i); %W=1.5cm/sec
    %transport_mod=102;
    [velocity_now,bed,moving]=transport103(velocity_now,bed,moving,i); %W=3cm/sec
    transport_mod=103;
    %[velocity_now,bed,moving]=transport104(velocity_now,bed,moving,i); %W=9cm/sec
    %transport_mod=104;
    %%%%%%%%%%%%%
    %Ribberink 1998
    %[velocity_now,bed,moving]=transport200(velocity_now,bed,moving,i,A,T,n,p);
    %transport_mod=200;
    %%%%%%%%%%%%%
    %Bailard transport 7/27/07
    %[velocity_now,bed,moving]=transport400(velocity_now, bed, moving, i, n, p, xs, ys);
    %transport_mod=400;
    
    %%%%%%%%%%%%
    %Call the subroutine that moves the sand that is in suspension
    [velocity_now, moving]=movesandinsus100(velocity_now, moving, i, jump_frac, n, p);
    
    %%%%%%%%%%%%%
    %Testing for angle of repose
    [bed, xs, ys]=angleofrepose100(bed,moving, n, p, yplot, xplot, i);
    
    %%%%%%%%%%%%%
    %Feedback
    [velocity_now,velocity_next, bed]=feedback101(velocity_now,velocity_next,bed, i, xs, ys, n, p);

    %%%%%%%%Other stuff%%%%%%%%
    %plot the new bed
    [bed]=plotbed(bed,velocity_now, moving, xplot, yplot, i); % to plot to the screen as the model runs (this expensive computationally)
    %F(i)=getframe; % to make a movie (this also expensive computationally)
    
    %calculate mass conservation and stats like rms and save to a file
    %[bed, moving, data]=bedinfo(bed,moving, A, S, T, var_amp, i, jump_frac,t,transport_mod);
    
    %calculate the spectrum of the bed and save to a file
    %[bed,spectra, plot_f]=bedspectra(bed, n, p, A, S, T, jump_frac, t, transport_mod,i);
    
    %This loop was used to save beds for making time series plots
    %write bed to a file after each 100 iterations
    %if cnt==100
    %    totcnt=totcnt+cnt;
    %    filename=['mripA75S20va100_bedt' num2str(totcnt)];
    %    save(filename, 'bed');
        %fid=fopen(['mripA75S20va100_bedt' num2str(totcnt)], 'w');
        %fwrite(fid,bed,'double');
        %fclose(fid);
    %    cnt=0;
    %else
    %    cnt=cnt+1;
    %end

end  %i loop through time


