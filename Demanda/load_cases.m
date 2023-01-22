 function [SLOAD,P,Q,carbus] = load_cases(mpc,a,hora,ncasos)
 
% By  Ing. Walter Villa
% Last modification Medellín 29.12.2015

% Purpose: this function create different scenarios of node load of power
% system, these scenarios are constructed based on three types of loads on
% the nodes of the power system and considering a normal distribution of the load on each node

% Inputs:
%   file => Filename of power system to evaluate
%   a => percentage deviation is applied to the mean (standard desviation)
%   hora => hour of admission to the demand curve
%   ncasos => number of load cases generated

% Outputs:
%    Sload: this corresponding to active and reactive power demand
%    P    : is a active power demand in each node according to demand curve to choose
%    Q    : is a reactive power demand in each node according to demand curve to choose
%    carbus: this variable is vector of the load buses.

% Calls: cargatipo

%   tipo => load type for each node in the system
%      1 => Residential demand
%      2 => Industrial demand
%      3 => Commercial demand

% Inicialización de la variable tipo
% tipo =[0 0 2 2 0 0 2 2 0 0 0 2 0 0 2 1 0 1 0 1 1 0 1 1 3 3 3 3 3 0 3 0 0 0 0 0 0 0 1];
% tipo =[0 0 1 1 0 0 1 1 0 0 0 1 0 0 1 1 0 1 0 1 1 0 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 0 1];
% tipo = [0 1 1 1 0 0 1 1 0 1 0 1 0 1 1 1 1 1 1 1 1 0 1 1 0 1 0 0 1 1]; % Sistema IEE 30 bus
% tipo = [0 1 1 1 1 1 0 0 1 1 1 1 1 1];  % Sistema IEE 14 bus

tipo = zeros(size(mpc.bus(:,1)));
carga_base = mpc.bus(:,3)+sqrt(-1)*mpc.bus(:,4); % S complex
carbus = find(mpc.bus(:,3)~=0);                  % PQ nodes
tipo(carbus,1)=1;
% carbus = find(tipo~=0);
SLOAD = zeros(ncasos,length(carbus)*2);
[P,Q] = cargatipo(tipo,hora,carga_base,carbus);
    
for flu = 1:ncasos
%     Num = floor((rand*(24)+1)); % otra forma Num = floor(8 + (24-8).*rand);
    Num = randi([1 24],1); % En este rango no se tiene en cuenta condiciones de baja demanda residencial
    Pcarga = random('normal',P(Num,:),P(Num,:)*a);
%     mpc.bus(carbus,3) = Pcarga';
    
    if all(Q(Num,:)>0)
        Qcarga = random('normal',Q(Num,:),Q(Num,:)*a);
    else
        neg = find(Q(Num,:)<0);
        Qcarga = random('normal',Q(Num,:),Q(Num,:)*a);
        Qcarga(neg) = random('normal',Q(Num,neg),Q(Num,neg)*-a);
    end
%     mpc.bus(carbus,4)= Qcarga';
    SLOAD(flu,:)=[Pcarga,Qcarga];
end

function [P,Q]= cargatipo(tipo,hora,PQ,carbus)
% Purpose: this function calculate active and reactive power according to the type of power demand in the nodes
% of power system. Types are Commercial,industrial and residential load
% Inputs: 
% Tipo:
%      1 => Residential demand
%      2 => Industrial demand
%      3 => Commercial demand
% Hora: corresponding to hour where calculated of conversion factor
% PQ: is the P+jQ load of node system.
% carbus: this variable is vector of the load buses.
% Output:
   % P : Active power for 24 hours
   % Q : Reactive power for 24 hours

% load demand curves of differents types of load in the system
Carga_Residencial= [69 65 62.1 56 58 61 64 76 90 95 98 100 99 99.5 105 97 96 94 93 92 91 90 88 72]'/100;
Carga_Industrial = [56 54 52 50 55 58 68 80 90 98 100 94 95 96 90 83 78 72 71 70 69 67 65 60 ]'/100;
Carga_Comercial  = [20 19 18.5 18 20 22 25 40 65 86 90 92 89 92 94 96 100 88 76 73 65 50 28 22]'/100;
% --------------------------------------------------------------------------
P = zeros(24,length(carbus));
Q = zeros(24,length(carbus));
fpico = 1.05;
Pd   = real(PQ(carbus));                    % Demandas nodales de potencia activa
Dsys = sum(Pd)*fpico;                       % Demanda total pico del sistema 
Pdm   = Pd*fpico;                           % Demandas nodales pico

for j = 1:length(carbus)   
switch tipo(carbus(j)) 
    case 1
        factor = Dsys/Carga_Residencial(hora);
        P(:,j) = Pd(j)*(1 + (factor*Carga_Residencial-Dsys)/Dsys)...  % Demanda inicial en cada hora para el nodo j
           .*(1 + Pdm(j)/factor*Carga_Residencial);
    case 2
        factor = Dsys/Carga_Industrial(hora);
        % perfil_Pd = Carga_Industrial*factor;  
        P(:,j) = Pd(j)*(1 + (Carga_Industrial*factor-Dsys)/Dsys)...
            .*(1 + Pdm(j)/factor*Carga_Industrial);
    case 3
        factor = Dsys/Carga_Comercial(hora);
        P(:,j) = Pd(j)*(1 + (factor*Carga_Comercial-Dsys)/Dsys)...
            .*(1 + Pdm(j)/factor*Carga_Comercial);
    otherwise
        disp('Unknown method')
end

if P(:,j)
    fp = real(PQ(carbus(j)))/abs(PQ(carbus(j)));
    if imag(PQ(carbus(j)))<0
        Q(:,j) = -P(:,j)*sqrt(1-fp^2)/fp;
    else
        Q(:,j) = P(:,j)*sqrt(1-fp^2)/fp;
    end
end
end
%     opt = mpoption('VERBOSE', 0, 'OUT_SYS_SUM', 0, 'OUT_BUS', 0, 'OUT_BRANCH', 0, 'OUT_ALL_LIM',0, 'OUT_V_LIM',0);
%     R = runopf(mpc,opt);
%     SGEN(flu,:)=[R.gen(1:end-G_ENS,2)',R.gen(1:end-G_ENS,3)'];
%     SLOAD(flu,:)=[R.bus(carbus,3)',R.bus(carbus,4)'];
%     Hora_opf(flu,:)= Num;
%     Pgenerada(:,flu)= R.gen(:,2);
%     BusGen = sort(R.gen(1:end-G_ENS,1));
%     %Vgen=zeros(length(BusGen),1);
%     Sgen = zeros(length(BusGen),1);
%     for hh=1:length(BusGen)
%         Loc=find(R.gen(:,1)==BusGen(hh));
%         Sgen(hh,1)=(R.gen(Loc,2)+R.gen(Loc,3)*sqrt(-1))/Sbase;
%     end
%     Vgen = R.bus(BusGen,8);                                  % Voltage magnitude of PV nodes
%     Agen = R.bus(BusGen,9)*pi/180;                           % Voltage angle of PV nodes
%     Vfasor = Vgen.*(cos(Agen)+sin(Agen)*sqrt(-1));           % Phasor of PV nodes
%     Vcar = R.bus(carbus,8);                                  % Voltage Magnitude of PQ nodes
%     Scar = (R.bus(carbus,3)+R.bus(carbus,4)*sqrt(-1))/Sbase; % Apparent Demand of PQ nodes
%     Zcarga(flu,:)= Vcar.^2./conj(Scar);
% end
% Pgmax = max(Pgenerada,[],2);
% Pmean = mean(Pgenerada,2);
%----------------------------------------------------------------------
    