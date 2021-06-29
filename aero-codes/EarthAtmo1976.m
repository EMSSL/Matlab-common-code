% Function to return the average temperature, pressure, density, and dynamic 
% viscosity for Earth's atmosphere for a given altitude input. The basis 
% for this model is the 1976 Standard Atmosphere Model. Written
% by Rajmohan Waghela for EMSSL. Modified by Christopher D. Yoder 
% on 05/21/2016.
%
% INPUT:
%   GEOALT      =   [m] the desired altitude between 0 and 105000 meters
% 
% OUTPUT:
%   T           =   [deg K] the temperature of the atmosphere
%   P           =   [Pa] the pressure of the atmosphere
%   RHO         =   [kg/m3] the density of the atmosphere
%   MU          =   [N-s/m2] the dynamic viscosity of the atmosphere

function [T,P,RHO,MU,R] = EarthAtmo1976(GEOALT)
% ATMO USES GEOALT TO CALCULATE THE ATMOSPHERIC PARAMETERS USING 
% 1976 U.S. STANDARD ATMOSHPERE 
% GEOALT - GEOMETRIC MSL ALTITUDE, meters
% T - TEMPERATURE, Kelvin
% P - PRESSURE, Pascal
% RHO - DENSITY, kg/m^3
% MU - VISCOSITY, N.sec/m^2
%   VISCOSITY IS CALCULATED USING SUTHERLAND'S EQUATION
% R - AIR SPECIFIC GAS CONSTANT, J/kg-K

% rev1 = CDY, 01-23-2021, added Rsp output 


% establish constants for the atmo model
% SEA LEVEL CONDITIONS
Ps = 1.01325*10^5;          % [Pa] 
RHOs = 1.2250;              % [kg/m^3]
Ts = 288.16;                % [deg K] 

% VISCOSITY PARAMETERS
lambda=1.512041288e-6;      % [Pa.sec.K^0.5]
C=120;                      % [deg K] 

% General Constants
g0=9.81;                    % [m/s^2]
R=287;                      % [N.m/(kg.K)]


% If STRUCTURE TO SELECT ZONE
if GEOALT<0
    % % ERROR MESSAGE 
    % error('Invalid Altitude Input. Altitude is below 0 meters.')
    warning('Altitude input is below 0 meters. Parameters for 0 meters are taken.')
    % GRADIENT REGION
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    
    % TEMPERATURE
    h=GEOALT;
    T=Ts+a1*(0);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
elseif GEOALT>=0 && GEOALT<11000 
    % GRADIENT REGION
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    
    % TEMPERATURE
    h=GEOALT;
    T=Ts+a1*(h);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
    
elseif GEOALT>=11000 && GEOALT<25000
    % ISOTHERMAL REGION
    
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    h=11000;
    % TEMPERATURE
    T=Ts+a1*(h);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    
    % Region 2 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=11000;
    h=GEOALT;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
    
elseif GEOALT>=25000 && GEOALT<47000
    % GRADIENT REGION
    
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    h=11000;
    % TEMPERATURE
    T=Ts+a1*(h);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    
    % Region 2 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=11000;
    h=25000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 3 Parameters
    a3=3*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=25000;
    h=GEOALT;
    Ts=T;
    T=Ts+a3*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a3*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a3*R)));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
    
elseif GEOALT>=47000 && GEOALT<53000
    % ISOTHERMAL REGION
    
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    h=11000;
    % TEMPERATURE
    T=Ts+a1*(h);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    
    % Region 2 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=11000;
    h=25000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 3 Parameters
    a3=3*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=25000;
    h=47000;
    Ts=T;
    T=Ts+a3*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a3*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a3*R)));
    
    % Region 4 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=47000;
    h=GEOALT;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
    
elseif GEOALT>=53000 && GEOALT<79000
    % GRADIENT REGION
    
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    h=11000;
    % TEMPERATURE
    T=Ts+a1*(h);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    
    % Region 2 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=11000;
    h=25000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 3 Parameters
    a3=3*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=25000;
    h=47000;
    Ts=T;
    T=Ts+a3*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a3*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a3*R)));
    
    % Region 4 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=47000;
    h=53000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 5 Parameter
    a5=-4.5*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=53000;
    h=GEOALT;
    Ts=T;
    T=Ts+a5*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a5*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a5*R)));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
    
elseif GEOALT>=79000 && GEOALT<90000
    % ISOTHERMAL REGION
    
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    h=11000;
    % TEMPERATURE
    T=Ts+a1*(h);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    
    % Region 2 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=11000;
    h=25000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 3 Parameters
    a3=3*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=25000;
    h=47000;
    Ts=T;
    T=Ts+a3*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a3*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a3*R)));
    
    % Region 4 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=47000;
    h=53000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 5 Parameter
    a5=-4.5*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=53000;
    h=79000;
    Ts=T;
    T=Ts+a5*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a5*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a5*R)));
    
    % Region 6 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=79000;
    h=GEOALT;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
    
elseif GEOALT>=90000 && GEOALT<=105000
    % GRADIENT REGION
    
    % Region 1 Parameters
    a1=-6.5*10^-3; % Kelvin/meter
    h=11000;
    % TEMPERATURE
    T=Ts+a1*(h);
    % DENSITY
    RHO=RHOs*(T/Ts)^(-(g0/(a1*R)+1));
    % PRESSURE
    P=Ps*(T/Ts)^(-(g0/(a1*R)));
    
    % Region 2 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=11000;
    h=25000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 3 Parameters
    a3=3*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=25000;
    h=47000;
    Ts=T;
    T=Ts+a3*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a3*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a3*R)));
    
    % Region 4 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=47000;
    h=53000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 5 Parameter
    a5=-4.5*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=53000;
    h=79000;
    Ts=T;
    T=Ts+a5*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a5*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a5*R)));
    
    % Region 6 Parameters
    % TEMPERATURE
    T=T;
    % DENSITY
    h1=79000;
    h=90000;
    RHO=RHO*exp(-(g0/(R*T))*(h-h1));
    % PRESSURE
    P=P*exp(-(g0/(R*T))*(h-h1));
    
    % Region 7 Parameter
    a7=4*10^-3; % Kelvin/meter
    % TEMPERATURE
    h1=90000;
    h=GEOALT;
    Ts=T;
    T=Ts+a7*(h-h1);
    % DENSITY
    RHO=RHO*(T/Ts)^(-(g0/(a7*R)+1));
    % PRESSURE
    P=P*(T/Ts)^(-(g0/(a7*R)));
    % VISCOSITY
    MU=lambda*T^(3/2)/(T+C);
    
elseif GEOALT>105000
    % BEYOND RANGE
    error('Invalid Altitude Input. Altitude is above 105000 meters.')
else
    error('DEFINED ERROR - Invalid Input - %d',GEOALT)
end

end