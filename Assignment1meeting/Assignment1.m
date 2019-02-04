%Assignment 1, Peter Al Ahmar, 100961570

clear;
close all;


q = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;           % Dirac constant
h = hb * 2 * pi;              % Planck constant
m = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;             % Boltzmann constant
eps = 8.854187817e-12;          % vacuum permittivity
mu = 1.2566370614e-6;           % vacuum permeability
c = 299792458;                  % speed of light
g = 9.80665;                    %metres (32.1740 ft) per s²
am = 1.66053892e-27;

iter = 1000;
T=300;
AtomSpacing = 0.5430710e-9;
nParticles = 50;
em= 0.26*m; %effective mass
vth = sqrt(kb*T/em);

strdDev = vth/sqrt(2);

xboundary = 200e-9;
yboundary = 200e-9;

initialx = rand(1,nParticles) .*xboundary;
initialy = rand(1,nParticles) .*yboundary;

velocityx = rand(1,nParticles) .*strdDev;
velocityy = rand(1,nParticles) .*strdDev;

dt = 10e-15;

temparray = zeros(1,iter); % temp array
temp = (1:1:iter);
figure(1);
plot(initialx,initialy,'.');
hold on;

initialS = zeros(1,nParticles);

tDiffer = 0;
SumtDiffer = 0;

col = 0;
pScatter = 1- (exp((-1*dt) / (0.2e-12))); %scatter 

for i=1: iter
    
    initialS1 = pScatter > rand(1,nParticles);
    velocityx(initialS1) = randn .* strdDev;
    velocityy(initialS1) = randn .* strdDev;
    
    if initialS(1) ~= initialS1(1)
       col=col +1;
       SumtDiffer = SumtDiffer + tDiffer;
       tDiffer = 0;
       
    else
        tDiffer = tDiffer + 1;
        
    end
    
    %x-boundary using periodic boundary condition
    xbn = initialx > 200e-9;
    xbn2 = initialx < 0;
    
    initialx(xbn) = initialx(xbn) - 200e-9;
    initialx(xbn2) = initialx(xbn2) + 200e-9;
    
    
    %y-boundary reflection
    
    ybn = initialy > 200e-9;
    ybn2 = initialy < 0;
    
    velocityy(ybn) = - velocityy(ybn);
    velocityy(ybn2) = - velocityy(ybn2);
    
    
    %update position 
    
    Nextx= initialx + velocityx .* dt;
    Nexty= initialy + velocityy .* dt;
    initialx = Nextx;
    initialy = Nexty;
    
    
    vrms = sqrt((velocityx .^ 2) + (velocityy .^ 2));
    temperature = (sqrt(2) * (mean(vrms)) ^2 *em) / kb; 
    temparray(1,i) = temperature;
    
    
    plot(initialx,initialy,'.');
    hold on;
    title(['Average Temperature: ', num2str(temperature)]);
    xlim([0,xboundary]);
    ylim([0,yboundary]);
    pause(0.1);
    
    hold on; 
    
    %ADDING box to the figure --- not working --
	% annotation('textbox',[yboundary,xboundary, yboundary/3, yboundary/3],'String','');
    
    
    line([0.75*xboundary/2 0.75*xboundary/2], [yboundary 2*yboundary/3]);
    line([1.05*xboundary/2 1.05*xboundary/2], [yboundary 2*yboundary/3]);
    line([0.75*xboundary/2 1.05*xboundary/2], [yboundary yboundary]);
    line([0.75*xboundary/2 1.05*xboundary/2], [2*yboundary/3 2*yboundary/3]);

    line([0.75*xboundary/2 0.75*xboundary/2], [0 yboundary/3]);
    line([1.05*xboundary/2 1.05*xboundary/2], [0 yboundary/3]);
    line([0.75*xboundary/2 1.05*xboundary/2], [0 0]);
    line([0.75*xboundary/2 1.05*xboundary/2], [yboundary/3 yboundary/3]);
    %Got help to create the boxes 
    %However I tried using the same logic that was used in PA4 to add
    %boundaries to them and only let the electrons flow outside the box but
    %was unable to do it correctly :/
    
end

%plot(X(nParticles),Y(nParticles));

%MFT is calculated below 
MeanFT=(SumtDiffer * dt)/col;

%Mean Free Path is calculated below
MeanFP=mean(vrms) *MeanFT;

figure(2);
histogram(vrms,10);
title('Thermal Velocity');

figure(3);
plot(temp,temparray);
title('Average Temperature')
xlabel('Time (s)');
ylabel('Temperature');

fprintf("Mean Free Path is %f", MeanFP);
fprintf("Mean Free Time is %f", MeanFT);
%ADDING box to the figure --- not working --



