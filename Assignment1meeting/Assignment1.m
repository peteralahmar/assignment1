%Assignment 1, Peter Al Ahmar, 100961570

clear;
close all;


q = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;           % Dirac constant
h = hb * 2 * pi;                % Planck constant
m = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;             % Boltzmann constant
eps = 8.854187817e-12;          % vacuum permittivity
mu = 1.2566370614e-6;           % vacuum permeability
c = 299792458;                  % speed of light
g = 9.80665;                    %metres (32.1740 ft) per s²
am = 1.66053892e-27;

iter = 200;
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

dt = 2* 10e-15;

temparray = zeros(1,iter); % temp array
temp = (1:1:iter);


tDiffer = 0;
SumtDiffer = 0;
col = 0;


inside_box = true;
while inside_box == true
    inside = ((initialx <= (1.15 * xboundary/2) & (initialx >= (0.85 * xboundary/2))) & ((initialy < (yboundary/3) | initialy >= (2*yboundary/3))));
    
    if (sum(inside) >0)
        initialx(inside) = rand(1,sum(inside)).* xboundary;
        initialy(inside) = rand(1,sum(inside)) .* yboundary;
    else
        inside_box = false;
    end
end

velocityx = randn(1,nParticles) .*strdDev;
velocityy = randn(1,nParticles) .*strdDev;


vrms = sqrt((velocityx .^ 2) + (velocityy .^ 2));
vrms_array =zeros(1,1000);

pScatter = 1- (exp((-1*dt) / (0.2e-12))); %scatter 
initialS = zeros(1,nParticles);


boundary = 0;

for i=1: iter
    
    initialS1 = pScatter > rand(1,nParticles);
    velocityx(initialS1) = randn(1,length(velocityx(initialS1))) .* strdDev;
    velocityy(initialS1) = randn(1,length(velocityy(initialS1))) .* strdDev;
    
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
    
    ybn = initialy >= 200e-9;
    ybn2 = initialy <= 0;
    
    velocityy(ybn) = - velocityy(ybn);
    velocityy(ybn2) = - velocityy(ybn2);
    
    oldx = initialx;
    oldy = initialy;
    
    
    inside = ((( initialx <(1.15 * xboundary/2)) & (initialx > (0.85 * xboundary/2))) & (( initialy< (yboundary/3)) | (initialy > (2*yboundary/3))));
    
    if (( boundary ==0) && (sum(inside) >=1))
        
        if (( oldx < (1.15 * xboundary/2)) & (oldx > (0.85 * (xboundary/2)) &(sum(inside) >=1)))
            
            if (initialy(inside) > (2* boundary/3))
                
                initialy(inside) = initialy(inside) - (2* (initialy(inside) - (2*yboundary/3)));
    
            elseif (initialy(inside) < ( boundary/3))
                
                initialy(inside) = initialy(inside) + (2* ((yboundary/3) - initialy(inside)));
                
            end
            
            velocityy(inside) = -velocityy(inside);
            initialx(inside) = initialx(inside) + (velocityx(inside) .*dt);
            initialy(inside) = initialy(inside) + (velocityy(inside) .*dt);
            
        else
            velocityx(inside) = -velocityx(inside);
            initialx(inside) = initialx(inside) + (velocityx(inside) .*dt);
            initialy(inside) = initialy(inside) + (velocityy(inside) .*dt);  
        end
        
        vrms = sqrt((velocityx .^ 2) + (velocityy .^ 2));
        
    elseif (( boundary ==1) && (sum(inside) >=1))
        
        if (( oldx < (1.15 * xboundary/2)) & (oldx > (0.85 * (xboundary/2)) & (velocityy(inside) > 0)))
          
            initialy(inside) = initialy(inside) - (2* (initialy(inside) - (2*yboundary/3)));
            velocityx(inside) = randn .* strdDev;
            velocityy(inside) = -1. * (abs(randn .* strdDev));
                
        elseif (( oldx < (1.15 * xboundary/2)) & (oldx > (0.85 * (xboundary/2)) & (velocityy(inside) < 0)))
            
            initialy(inside) = initialy(inside) - (2* ((boundaryy/3) - initialy(inside)));
            velocityx(inside) = randn .* strdDev;
            velocityy(inside) =(abs(randn .* strdDev));
            
        elseif (velocityx(inside) >0)
            
            initialx(inside) = initialx(inside) - (2* (initialx(inside) - (0.85*yboundarx/2)));
            velocityy(inside) = randn .* strdDev;
            velocityx(inside) = -1. * (abs(randn .* strdDev));
            
        else
                        
            initialx(inside) = initialx(inside) + ((2* (1.15*boundaryx/2) - initialx(inside)));
            velocityy(inside) = (abs(randn .* strdDev));
            velocityx(inside) = (abs(randn .* strdDev));
            
        end
        
        vrms = sqrt((velocityx .^ 2) + (velocityy .^ 2));
        
        
    end

    %update position 
    
    initialx= oldx + (velocityx .* dt);
    initialy= oldy + (velocityy .* dt);

    vrms = sqrt((velocityx .^ 2) + (velocityy .^ 2));
    
    

    temperature = (sqrt(2) * (mean(vrms)) ^2 *em) / kb; 
    temparray(1,i) = temperature;
    
    figure(1);
    plot(initialx,initialy,'.');
    xlabel("x-Position");
    ylabel("y-Position");

    title(['Average Temperature: ', num2str(temperature)]);

    xlim([0,xboundary]);
    ylim([0,yboundary]);
    hold on; 

    
    line([0.85*xboundary/2 0.85*xboundary/2], [yboundary 2*yboundary/3]);
    line([1.15*xboundary/2 1.15*xboundary/2], [yboundary 2*yboundary/3]);
    line([0.85*xboundary/2 1.15*xboundary/2], [yboundary yboundary]);
    line([0.85*xboundary/2 1.15*xboundary/2], [2*yboundary/3 2*yboundary/3]);

    line([0.85*xboundary/2 0.85*xboundary/2], [0 yboundary/3]);
    line([1.15*xboundary/2 1.15*xboundary/2], [0 yboundary/3]);
    line([0.85*xboundary/2 1.15*xboundary/2], [0 0]);
    line([0.85*xboundary/2 1.15*xboundary/2], [yboundary/3 yboundary/3]);
    
end

%plot of Temperature vs time

figure(2);
plot(temp,temparray);
title('Average Temperature')
xlabel('Time (s)');
ylabel('Temperature');
hold on


%MFT is calculated below 
MeanFT=(SumtDiffer * dt)/col;

%Mean Free Path is calculated below
MeanFP= mean(vrms) *MeanFT;

fprintf("Mean Free Path is %f", MeanFP);
fprintf("Mean Free Time is %f", MeanFT);

[xgr,ygr] = meshgrid(0:(xboundary/10):xboundary, 0:(yboundary/10):yboundary);
electron_M = zeros(11,11);
temperature_M = zeros(11,11);
num_Elec = 0;
Total_vel = 0;

for k = 1:10
    x_min = xgr(1,k);
    x_max = xgr(1,k+1);
    
    for m = 1:10    
        y_min = ygr(m,1);
        y_max = ygr(m+1,1);
        
        for n = 1:nParticles
            if ((initialx(n) >x_min) && (initialx(n) < x_max) && ((initialy(n) > y_min) && initialy(n) <y_max))
                num_Elec = num_Elec +1;
                Total_vel = Total_vel + sqrt((velocityx(n) .^2) + (velocityy(n) .^2));
                electron_M(k,m) = electron_M(k,m) + 1;
                temperature_M(k,m) = ((sqrt(2)*(Total_vel/num_Elec)^2)*em)/kb;
            end
             
        end
        Total_vel = 0;
        num_Elec = 0;
        
    end
    
end

%Thermal Velocity Histogram
figure(3);
histogram(vrms,10);
title('thermal Velocity Histogram');
    
% Creating an Electron Map
figure(4);
surf(electron_M);
title('Electron Density Map');
    
%Creating a Temperature Map
figure(5);
surf(temperature_M);
title('Temperature Map');



