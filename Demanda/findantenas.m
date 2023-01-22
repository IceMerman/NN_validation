function [busidx,lineidx] = findantenas(mpc)

% This function realize topological analysis of power system, it
% identifies buses that are connected to the network through a single
% line that can causes no convergence of the power flow algorithm.


% FINDANTENAS finds buses in antenna and returns the indexes
%              of buses and of the unique connected lines
% Inputs:  mpc  power system matrix structure
% BUSIDX:  indexes of buses in antena
% LINEIDX: indexes of lines connected to buses in antena
%


busidx = [];
lineidx = 0;

% ---------------------------
nb = size(mpc.bus,1);
n = size(mpc.branch,1);
nl = [1:n]';

fr = mpc.branch(:,1);                     % Numeración externa (sistema original)
to = mpc.branch(:,2);                     % Numeración externa (sistema original)

% -------- Revisión de la numeración de interna de la matrices -------------
% if max(bus.con(:,1))> nb
    % pre allocation
%     Fr = zeros(size(fr,1),1);
%     To = zeros(size(to,1),1);
%     idx_b =  bus.con(:,end);        % numeración interna de los nodos
%     dentro de la matrix bus
%     for i=1:idx_b(end)
%         Ir = find(fr == bus.con(i,1));
%         It = find(to == bus.con(i,1));
%         Fr(Ir)= i;
%         To(It)= i;
%     end
% else
%     Fr = fr;
%     To = to;
% end

% Otra forma de hacer la identificación de los nodos con la numeración interna
busmax = max(mpc.bus(:,1));
bus_int = zeros(busmax,1);
ibus = (1:nb)';
bus_int(round(mpc.bus(:,1))) = ibus;
Fr = bus_int(fr);
To = bus_int(to);

ivec = sparse(nl,Fr,1,n,nb);
jvec = sparse(nl,To,1,n,nb);
% A = ivec - jvec;                       % Matriz de Incidencia o de conexiones
busidx = (find(sum([ivec;jvec],1)==1))';
%  % sparse matrix formulation
%   tap = ones(n,1);
%   iline = [1:1:nline]';
%   C_from = sparse(Fr,nl,tap,nb,n,n);
%   C_to = sparse(To,nl,ones(n,1),nb,n,n);
%   C_line = C_from - C_to;

if isempty(busidx)
    lineidx = 0;
%     disp('All lines are used for (N-1) contingency evaluations.')
else
%     lineidx = find(sum(ivec(:,busidx)+jvec(:,busidx),2));
    [lineidx,ff] = find(ivec(:,busidx)+jvec(:,busidx));
end
% -------------------------------------------------------------
% se asegura de no sacar el generador Slack del sistema
%  slack = find(bus.con(busidx,2) == 3);   
%  busidx(slack,:)= [];
%  lineidx(slack,:) = [];
 % ------------------------------------------------------------
 