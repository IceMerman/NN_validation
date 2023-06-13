%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%         DATABASE FOR IA USING MCS BASED OF PROBABILISTIC POWER FLOW
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% By Ing. Walter Villa
% San Juan, 23.01.2013

% Ultima modificacion:
% Medellin, 31.10.2022

% Comentarios:
% No se tiene en cuenta el modelado de Generadores de Energia no Suministrada.
% Se incluye la salida de generadores, transformadores y lineas de transmision.
% Se considera el analisis de nodos Antenna (PV y PQ).
% Se realiza analisis de contingencia
%%

% close all
clear all
clear classes
clc
% tic

%% -------------------------------------------------------------------------
%**************************************************************************
%                               Reading PDFs
%**************************************************************************
N = 100; % Amount of random numbers to be generated

compute_contigencies = 1;

%% Load data of power system (Matpower Format)
mpc = NewEngland_New;
file = 'New_England';

%% Identification of load nodes
nonnull_pload_idx = find(mpc.bus(:,3)~=0); % identificanco los nodos de carga PQ en sistema bajo prueba

%% Load demand cases (generacion de escenarios de demanda de forma probabilistica)

% SLOAD = load_cases(mpc,a,hora,N)
% SLOAD: corresponde a potencia activa y reactiva calculada para cada nodo PQ
% mpc: nombre del sistema a evaluar
% a: Desviacion de demanda considerada para distribucion normal de la demanda en cada nodo
% hora: Momento de admision a la curva de demanda
% N: numero de casos construidos (escenarios de demanda)
%**************************************************************************
SLOAD = load_cases(mpc,0.08,9,N);    %24 hours
%%**************************************************************************

%% -------------------------------------------------------------------------
%                 (N-1) Contigencies selection by power system
% -------------------------------------------------------------------------
%  Topological analysis finding "antenas": se identifican nodos
%  antenas(nodos que se conecta la sistema atraves de una rama ya sea LT o
%  transformador).
[antenas,fromto] = findantenas(mpc);

% Transmision lines:
Tlines_idx = find(mpc.branch(:,9)==0); % numeracion interna de la matriz branch

% Transformers:
Transf_idx = find(mpc.branch(:,9)~= 0); % numeracion interna de la matriz branch

% Seleccion de lineas transmision y transformadores del sistema de potencia

list_LT = sort([Tlines_idx;Transf_idx]) ;
nLine = size(list_LT);                           % # de ramas del sistema que pueden salir en la contigencia (Lines and transformers)

% Forced outage rate (q) of transmission lines, q = 0.04
qtl = (0.04);            % ones(1,length(Tlines_idx))*(0.2/100);
P0 = 1-exp(-1*qtl);      % probabilidad de contingencia N-1 en el sistema. Tomado del Libro: Probabilistic Transmission system planning

poperationL = random('Uniform',0,1,[N length(list_LT)]);  % Uniform distribution for transmission lines outages
% =======================================================================================================================
% Generators:
% Seleccion de generadores del sistema de potencia (Seleccionar unicamente  tipo PV)
fgen = find(mpc.bus(:,2) == 2);                 % generadores (PV) del sistema #s internos de la matriz
ngen = size(fgen,1);                            % # de Generadores que se consideran en las contigencias. No el nodo Slack.
poperationG = random('Uniform',0,1,[N ngen]);   % Uniform distribution for generator outages
nconti = nLine + ngen;                          % # de contingencias que se puede tener en cuenta

%%
%**************************************************************************
%               Monte Carlo-based Probabilistic optimal load Flow
%**************************************************************************

% numero de corridas del Monte Carlo
numeval = N;
% preallocation de las base de datos (estructura):
datosLoad = struct('loading',zeros(numeval,1));
OPF.VM   = zeros(size(mpc.bus,1),numeval);
OPF.VA   = zeros(size(mpc.bus,1),numeval);
OPF.Qgen = zeros(size(mpc.gen,1),numeval);
OPF.Pgen = zeros(size(mpc.gen,1),numeval);
OPF.PgenDC = zeros(size(mpc.gen,1),numeval);
% OPF.fl_pfrom = zeros(size(mpc.branch,1),numeval);
% OPF.fl_qfrom = zeros(size(mpc.branch,1),numeval);
% OPF.fl_pto = zeros(size(mpc.branch,1),numeval);
% OPF.fl_qto = zeros(size(mpc.branch,1),numeval);
issecure = zeros(1, N);
contigency = zeros(1, N);
branchOut = zeros(1, N);
genOut = zeros(1, N);

% Divergente cases:
% inizalization of counter
Ncd = 0;
I=0;
% Number of generators that acts as load shedding (must be that last generators in the mpc file)
n_shedding = 3;

% MATPOWER options
% opt = mpoption('VERBOSE', 0, 'OUT_SYS_SUM', 0, 'OUT_AREA_SUM', 0,'OUT_BUS', 0,'OUT_BRANCH', 0,'OUT_GEN', 0,'OUT_ALL_LIM',0);
%opt = mpoption('pf.enforce_q_lims', 0, 'opf.violation', 1e-4, 'opf.flow_lim', 'S', 'verbose', 0, 'out.all', 0);


iss = 0;                   % Inicialization of MCS counter
branch = mpc.branch;       % Initial configuration of Network (Base case)
gen = mpc.gen;             % Initial configuration of Generators (Base case)
bus = mpc.bus;             % Initial configuration of Buses (Base case)
totgen = sum(mpc.gen(:,2));

%%
while (iss <= numeval-1)
    iss = iss+1;         % update of counter
    % issecure(iss) = 0;   % Set the case secure flag as 0
    fprintf(' Run Monte Carlo simulation # %5d \n',iss);
    mpc.branch = branch;
    mpc.gen = gen;
    mpc.bus = bus;
    
    % Set line power limits to infinite (to solve a DC single node opf)
    %mpc.branch(:,6:8) = 1e20;
    % Randomize gen cost
    [rsize, csize] = size(mpc.gencost);
    mpc.gencost(1:end-n_shedding, 5) = rand(rsize-n_shedding,1);
    mpc.gencost(1:end-n_shedding, 6) = 100*rand(rsize-n_shedding,1);
    mpc.gencost(1:end-n_shedding, 7) = 15000*rand(rsize-n_shedding,1);
    
    %To update the generated random values of demand
    mpc.bus(nonnull_pload_idx,3) = SLOAD(iss,1:end/2)';      % Load active power
    mpc.bus(nonnull_pload_idx,4) = SLOAD(iss,end/2+1:end)';  % Load reactive power
    %%
    % Economic dispatch (DC optimal power flow) & BASE CASE
    opt = mpoption( 'opf.violation', 1e-4,'verbose', 0, 'out.all', 0);
    mpc.branch(:, 6) = 0; % disable line flow limits (mimic no network case) %-----------
    r = rundcopf(mpc, opt);
    %%
    % Revision de la convergencia del OPF
    % Get DC result for generators
    OPF.PgenDC(:,iss) = r.gen(:,2);
    mpc.gen(:,2) = r.gen(:,2);
    
    % Restore line limits
    %mpc.branch(:, 6:8) = branch(:, 6:8);
    mpc.branch(:, 6) = branch(:, 6);
    
    %Se debe cambiar la options para considerar los limites de potencia en
    % en generadores para el AC PowerFlow
    %---------------------------------------------------------------------
    opt = mpoption('pf.enforce_q_lims',1,'pf.alg','NR','verbose', 0, 'out.all', 0);
    %---------------------------------------------------------------------
    if r.success==0
        disp('The Economical Dispatch (DCOPF) diverges ');
        Ncd = Ncd + 1;
        issecure(iss) = 64;
    else
        check_status = checkACLimits(mpc, n_shedding, opt);
        issecure(iss) = issecure(iss) + check_status;
        if check_status>1
            Ncd = Ncd + 1;
        end
    end
    
    % N-1 CONTINGENCY SELECTION AND ANALYSIS
    % To update of the generators and lines in failure operation
    
    Il = ((poperationL (iss,:) <= P0))';
    Ig = ((poperationG (iss,:) <= P0))';
    
    if (logical(nnz(Ig))) && (logical(nnz(Il)))
        if logical(max(poperationL(iss,Il)) > max(poperationG(iss,Ig)))
            % Making sure to N-1 contingency
            I = find(poperationL (iss,:) == max(poperationL(iss,Il)));
            Ig =[];
        else
            % Making sure to N-1 contingency
            I = find(poperationG (iss,:) == max(poperationG(iss,Ig)));
            Il = [];
        end
        
    elseif logical(nnz(Ig))
        I = find(poperationG (iss,:) == max(poperationG(iss,Ig)));
        Il = [];
        
    elseif logical(nnz(Il))
        I = find(poperationL (iss,:) == max(poperationL(iss,Il)));
        Ig =[];
    end
    
    if isempty(Ig)
        if find(I, 1, 'first') == 22    % outages of this branch yields aislands
            I(I)= 0;
            if ~isempty(dfm), I(dfm(1))= 1;end
        end
        
        % change to status failure (outages TL)
        mpc.branch(I,11)= 0;
        branchOut(iss) = I;
        
        % Para las lineas que unen los nodos antena tenemos que analizar el tipo de nodo
        outline = find(I);
        if ~isempty(outline)
            idx = find(fromto == outline);
            if ~isempty(idx)
                tipo = mpc.bus(antenas(idx),2);     % tipo de nodo (# internos)
                if tipo == 3, continue; end;
                
                switch tipo                                 % se remueve el nodo de la antena (Para evitar nodos en isla)
                    case 1 % PQ node
                        mpc.bus(antenas(idx),:) = [];
                    case 2 % PV node
                        mpc.bus(antenas(idx),:) = [];
                        %                         ixdgen = mpc.gen(:,1);                                       % Numeracion original del sistema
                        fg = find(r.gen(:,1) == mpc.bus(antenas(idx),1));            % Antenas(idx));
                        Pgs = mpc.gen(fg,2);
                        mpc.gen(fg,:)= [];
                        if (Pgs~=0)
                            mpc.gen(:,2) = r.gen(:,2).*(1 + Pgs./totgen);            % distribuye la potencia que sale
                        end
                end
            end
        end
    end
    
    if isempty(Il)
        fg = find(mpc.gen(:,1) == mpc.bus(fgen(I),1),1,'first');
        % change to status failure (outages of generator)
        mpc.gen(fg,8)= 0;
        genOut(iss) = fg;
    end
    
    % N-1 Contingency analysis:
    if nnz(I) >= 1 && issecure(iss) == 1 && compute_contigencies == 1
%         r.gen(:, 2)= r.gen(:,2);
        %opt = mpoption('pf.enforce_q_lims', 1, 'verbose', 0, 'out.all', 0);
        qr = runpf(mpc,opt);
        contigency(iss)=1;
        
        % Revision de la convergencia del PF
        if qr.success == 0
            disp('The N-1 Power Flow diverges ');
            Ncd = Ncd + 1;
            issecure(iss) = issecure(iss) + 128;
        else 
            check_status = checkACLimits(qr, n_shedding, opt);
            if check_status > 1
                Ncd = Ncd + 1;
                issecure(iss) = issecure(iss) + check_status - 1;
            else
                r = qr;
            end
        end
        
    else
        %opt = mpoption('pf.enforce_q_lims', 1, 'verbose', 0, 'out.all', 0);
        r = runpf(mpc, opt);
    end
    %     [busout, genout, branchout, f, success, info, et, g, jac, xr, pimul] = opf(baseMVA, bus, gen, branch, areas, gencost,opt);
    
    %% To save voltages, angles and other importants variables in each case
    OPF.VM(:,iss)   = r.bus(:,8);
    OPF.VA(:,iss)   = r.bus(:,9);
    OPF.Qgen(:,iss) = r.gen(:,3);
    OPF.Pgen(:,iss) = r.gen(:,2);
    OPF.Qden(:,iss) = r.bus(nonnull_pload_idx,4);
    OPF.Pden(:,iss) = r.bus(nonnull_pload_idx,3);
    OPF.fl_pfrom(:,iss)= r.branch(:,14);
    OPF.fl_qfrom(:,iss)= r.branch(:,15);
    OPF.fl_pto(:,iss)= r.branch(:,16);
    OPF.fl_qto(:,iss)= r.branch(:,17);
    
end
% toc
%% Export data

fprintf('>> Safe ration: %.2f\n', sum(issecure>1)/length(issecure));

% namebase = ['BaseDatosNE_',file,'31_10_2022.mat'];
% namebase = ['BaseDatosNE_',file,datestr(today(), 'yyyymmdd', 'local'),'.mat'];
% save(namebase, 'OPF','poperationL','poperationG','SLOAD');

VM = OPF.VM;
VA = OPF.VA;
Qgen = OPF.Qgen;
Pgen = OPF.Pgen;
Qload = OPF.Qden;
Pload = OPF.Pden;
PgenDC = OPF.PgenDC;

dbname = sprintf('db%d',N);
save(dbname, 'VM', 'VA', 'SLOAD', 'issecure', 'Qgen', 'Pgen', 'contigency', 'branchOut', 'genOut', 'PgenDC');