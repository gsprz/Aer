close all
clear all
clc
 
% The script plots a graphic of C_L, C_D  in function of alpha and C_L vs
% C_D for the configuration with the wing only.
% If you need to add the tail, you can add the parameters at the end of the
% script and change config.NCorpi = 2.
%% Test Case 1
% We can start with the Cessna's graphics

U_Inf_Mag = 60;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;

config.NCorpi = 1;

config.RootChord = [1.625];
config.DihedralAngle = [1.44]; % [°]
config.SweepAngle = [0]; % [°]
config.TaperRatio = [0.672]; 
config.Span = [11];
config.AspectRatio = [7.32]; 
config.LEPosition_X = [0];
config.LEPosition_Y = [0];
config.LEPosition_Z = [0];

config.RotationAngle_X = [0];
% config.RotationAngle_Y = [5,-1.5];
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [5];
config.ChordwiseDiscr = [5];

% Impose of alpha variation
RotationY = linspace(-5,12,18);
for angleIndex = 1:length(RotationY)
   for iCorpo = 1:config.NCorpi
       config.RotationAngle_Y(iCorpo) = RotationY(angleIndex);
    end
%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2;
% % Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create the geometry structure

ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);

for iCorpo = 1:config.NCorpi

    % [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
    % [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure_rec(config, iCorpo);
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
end
    

%% Matrices initialization

NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matriceA = zeros(NPanelsTot, NPanelsTot);
TermineNoto = zeros(NPanelsTot, 1);

%% Construction of the matrix

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
   
            columnIndex = 0;
            
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            
            for jCorpo = 1:config.NCorpi
                
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        
                        % Update column index
                        columnIndex = columnIndex + 1;
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere);
                       
                        
                    end
                end

            end
            
        
            
        end
    end
end

%% Costruzione del termine noto
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
  
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);

            
        end
    end
end

%% Solve the linear system

Solution = linsolve(matriceA, TermineNoto);

Gamma = cell(config.NCorpi, 1);

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    Gamma{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
    
     % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
            
            Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
        end
        
    end
    
end


%% Compute the 2D and 3D Lift

% We can calculate the lift with Kutta - Joukowsky theorem.
% Calculation of 2D lift (lift of the panels)
L_2D = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
L_2D{iCorpo} = zeros(2*config.SemiSpanwiseDiscr(iCorpo), 1);
      for k = 1:2*config.SemiSpanwiseDiscr(iCorpo)
          for i = 1:config.ChordwiseDiscr(iCorpo)
            L_2D{iCorpo}(k) = L_2D{iCorpo}(k) + rho * U_Inf_Mag * Gamma{iCorpo}(i,k) * cos(deg2rad(config.DihedralAngle(iCorpo)));
          end
      end

disp(L_2D{iCorpo});
end
C_l = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
for k = 1:2 * config.SemiSpanwiseDiscr
C_l{iCorpo}(k) = L_2D{iCorpo}(k)/(0.5 * rho * U_Inf_Mag^2); % ogni profilo
end
end
% Calcoliamo la sezione trasversale dei singoli pannelli
Deltab = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
Deltab{iCorpo} = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
end
% Totale lift
L_compl = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
    L_compl{iCorpo} = 0;
    for  i = 1:2 * config.SemiSpanwiseDiscr(iCorpo)
         L_compl{iCorpo} = L_compl{iCorpo} + Deltab{iCorpo} * L_2D{iCorpo}(i);
    end
end
% Coefficient of lift
C_L = zeros(1,config.NCorpi);
for iCorpo = 1:config.NCorpi
C_L(iCorpo) = L_compl{iCorpo}/(0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));
C_L_totale_per_angle(angleIndex) = sum(C_L);
end
%% Alpha induced and Cd
% Preallocate induced velocity and angle of attack cells
v_induced = cell(config.NCorpi, 1); 
alpha_induced = cell(config.NCorpi, 1);

for iCorpo = 1:config.NCorpi
    % Preallocate U_total for this body
    idx = 1; % Initialize index counter for U_total

    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        for SpanPanel_i = 1:2 * config.SemiSpanwiseDiscr(iCorpo)

            InducedPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            U_total = [];

            for jCorpo = 1:config.NCorpi
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    for SpanPanel_j = 1:2 * config.SemiSpanwiseDiscr(jCorpo)

                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U_total_v = induced_velocity(Gamma{jCorpo}(ChordPanel_j, SpanPanel_j), InducedPointHere, Extreme_1, Extreme_2);

                        % Store induced velocity in preallocated array
                        U_total(:, idx) = sum(U_total_v, 2); 
                        idx = idx + 1;

                    end
                end
            end


            U = sum(U_total, 2); % Sum across columns
            v_normal = -dot(U, NormalHere); % Normal velocity
            v_induced{iCorpo}(ChordPanel_i, SpanPanel_i) = v_normal; 
            alpha_induced{iCorpo}(ChordPanel_i, SpanPanel_i) = -atand(v_normal / U_Inf_Mag);
        end
    end
end

% Initialize variables for results
D_2D = cell(config.NCorpi, 1); % Drag distribution for each body
C_d = cell(config.NCorpi, 1); % Drag coefficient distribution for each body
C_D = zeros(config.NCorpi, 1); % Total drag coefficient for each body

% Loop through each body
for iCorpo = 1:config.NCorpi
    % Initialize drag distribution for the current body
    D_2D{iCorpo} = zeros(2 * config.SemiSpanwiseDiscr(iCorpo), 1);

    % Loop through each spanwise segment
    for k = 1:2 * config.SemiSpanwiseDiscr(iCorpo)
        % Sum contributions across chordwise panels
        for i = 1:config.ChordwiseDiscr(iCorpo)
            D_2D{iCorpo}(k) = D_2D{iCorpo}(k) + ...
                rho * U_Inf_Mag * Gamma{iCorpo}(i, k) * ...
                cos(deg2rad(config.DihedralAngle(iCorpo))) * ...
                sind(-alpha_induced{iCorpo}(i, k));
        end

        % Compute the drag coefficient for the spanwise section
        C_d{iCorpo}(k) = D_2D{iCorpo}(k) / (0.5 * rho * U_Inf_Mag^2);
    end
end

D_compl = zeros(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
    % Compute total drag for the current body
    Deltab = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo); % Panel spanwise width

    for i = 1:2 * config.SemiSpanwiseDiscr(iCorpo)
        D_compl(iCorpo) = D_compl(iCorpo) + Deltab * D_2D{iCorpo}(i);
    end

    % Compute total drag coefficient for the current body
    C_D(iCorpo) = D_compl(iCorpo) / (0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));
end
    C_D_totale_per_angle(angleIndex) = sum(C_D);
end

% Plot of CL
figure(1);
hold on;
plot(RotationY,C_L_totale_per_angle,'Color',[0.7 0 0],'LineWidth',2);
xlabel('\alpha (°)', 'FontSize',20,'FontWeight','bold');
ylabel('Coefficiente di Portanza','FontSize',20,'FontWeight','bold');
grid on;
set(gca,'FontSize',20,'FontWeight','normal');
hold off;

% Plot of CD
figure(2);
hold on;
plot(RotationY,C_D_totale_per_angle,'Color',[0 0 0.6],'LineWidth',2);
xlabel('\alpha (°)','FontSize',20,'FontWeight','bold');
ylabel('Coefficiente di Resistenza', 'FontSize',20,'FontWeight','bold');
grid on;
set(gca,'FontSize',20,'FontWeight','normal');
hold off

% Plot of CD vs CL
figure(3);
hold on;
plot(C_D_totale_per_angle,C_L_totale_per_angle,'Color',[0.2,0.6,0.2],'LineWidth',2);
xlabel('Coefficiente di Resistenza, C_D', 'FontSize',20,'FontWeight','bold');
ylabel('Coefficiente di Portanza, C_L','FontSize',20,'FontWeight','bold');
grid on;
set(gca,'FontSize',20,'FontWeight','bold');
hold off

% Same thing for Diamond DA40 to compare results
U_Inf_Mag = 60;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;

config.NCorpi = 1;

config.RootChord = [1.53];
config.DihedralAngle = [5]; % [°]
config.SweepAngle = [1]; % [°]
config.TaperRatio = [0.647]; 
config.Span = [11.94];
config.AspectRatio = [10.53]; 
config.LEPosition_X = [0];
config.LEPosition_Y = [0];
config.LEPosition_Z = [0];

config.RotationAngle_X = [0];
% config.RotationAngle_Y = [5,-1.5];
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [5];
config.ChordwiseDiscr = [5];

RotationY = linspace(-5,12,18);
for angleIndex = 1:length(RotationY)
   for iCorpo = 1:config.NCorpi
       config.RotationAngle_Y(iCorpo) = RotationY(angleIndex);
    end
%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2;
% % Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create the geometry structure

ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);

for iCorpo = 1:config.NCorpi

    % [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
    % [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure_rec(config, iCorpo);
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
end
    

%% Matrices initialization

NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matriceA = zeros(NPanelsTot, NPanelsTot);
TermineNoto = zeros(NPanelsTot, 1);

%% Construction of the matrix

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
   
            columnIndex = 0;
            
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            
            for jCorpo = 1:config.NCorpi
                
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        
                        % Update column index
                        columnIndex = columnIndex + 1;
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere);
                       
                        
                    end
                end

            end
            
        
            
        end
    end
end

%% Costruzione del termine noto
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
  
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);

            
        end
    end
end

%% Solve the linear system

Solution = linsolve(matriceA, TermineNoto);

Gamma = cell(config.NCorpi, 1);

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    Gamma{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
    
     % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
            
            Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
        end
        
    end
    
end

%nPanelsCorpo1 = 2 * config.SemiSpanwiseDiscr(1);
%nPanelsCorpo2 = 2 * config.SemiSpanwiseDiscr(2);

%span_corpo1 = linspace(-config.SemiSpan(1), config.SemiSpan(1), nPanelsCorpo1);
%span_corpo2 = linspace(-config.SemiSpan(2), config.SemiSpan(2), nPanelsCorpo2);

%Gamma_totale_corpo1 = Gamma{1};
%Gamma_totale_corpo2 = Gamma{2};

%Gamma_interpolata_corpo2 = zeros(config.ChordwiseDiscr(1), 2 * config.SemiSpanwiseDiscr(1));
%for k = 1:config.ChordwiseDiscr(1)
    %for j = 1:2 * config.SemiSpanwiseDiscr(1)
        %if span_corpo1(j) >= -config.SemiSpan(2) && span_corpo1(j) <= config.SemiSpan(2)
    %Gamma_interpolata_corpo2(k,j) = spline(span_corpo2, Gamma_totale_corpo2(k,:), span_corpo1(j));
        %else
            %Gamma_interpolata_corpo2(k,j) = 0;
        %end
    %end
%end
%Gamma_totale_combined = Gamma_totale_corpo1 + Gamma_interpolata_corpo2;
%Gamma_totale_combined_somma = sum(Gamma_totale_combined,1);
% Grafico della distribuzione di circolazione lungo lo span
%figure;
%plot(span_corpo1, Gamma_totale_combined_somma, '-', 'LineWidth',2, 'Color',[0.9,0.7,0]);
%xlabel('Posizione lungo lo span, y','FontSize',12,'FontWeight','bold');
%ylabel('Circolazione, \Gamma','FontSize',12,'FontWeight','bold');
%title('Distribuzione di Circolazione Totale lungo lo Span','FontSize',14,'FontWeight','bold');
%grid on;
%set(gca,'FontSize',10,'FontWeight','normal');
%hold off



%% Compute the 2D and 3D Lift

% Possiamo calcolare il lift attraverso il teorema di Kutta - Joukowski.
% Calcoliamo per prima cosa il contributo di portanza 2D, ovvero dei
% singoli pannelli in cui abbiamo suddiviso l'apertura alare.
L_2D = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
L_2D{iCorpo} = zeros(2*config.SemiSpanwiseDiscr(iCorpo), 1);
      for k = 1:2*config.SemiSpanwiseDiscr(iCorpo)
          for i = 1:config.ChordwiseDiscr(iCorpo)
            L_2D{iCorpo}(k) = L_2D{iCorpo}(k) + rho * U_Inf_Mag * Gamma{iCorpo}(i,k) * cos(deg2rad(config.DihedralAngle(iCorpo)));
          end
      end

disp(L_2D{iCorpo});
end
C_l = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
for k = 1:2 * config.SemiSpanwiseDiscr
C_l{iCorpo}(k) = L_2D{iCorpo}(k)/(0.5 * rho * U_Inf_Mag^2); % ogni profilo
end
end
% Calcoliamo la sezione trasversale dei singoli pannelli
Deltab = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
Deltab{iCorpo} = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
end
% Contributo di portata complessiva
L_compl = cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
    L_compl{iCorpo} = 0;
    for  i = 1:2 * config.SemiSpanwiseDiscr(iCorpo)
         L_compl{iCorpo} = L_compl{iCorpo} + Deltab{iCorpo} * L_2D{iCorpo}(i);
    end
end
C_L = zeros(1,config.NCorpi);
for iCorpo = 1:config.NCorpi
C_L(iCorpo) = L_compl{iCorpo}/(0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));
C_L_totale_per_angle(angleIndex) = sum(C_L);
end
%% Alpha induced and Cd
% Preallocate induced velocity and angle of attack cells
v_induced = cell(config.NCorpi, 1); 
alpha_induced = cell(config.NCorpi, 1);

for iCorpo = 1:config.NCorpi
    % Preallocate U_total for this body
    idx = 1; % Initialize index counter for U_total

    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        for SpanPanel_i = 1:2 * config.SemiSpanwiseDiscr(iCorpo)

            InducedPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            U_total = [];

            for jCorpo = 1:config.NCorpi
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    for SpanPanel_j = 1:2 * config.SemiSpanwiseDiscr(jCorpo)

                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U_total_v = induced_velocity(Gamma{jCorpo}(ChordPanel_j, SpanPanel_j), InducedPointHere, Extreme_1, Extreme_2);

                        % Store induced velocity in preallocated array
                        U_total(:, idx) = sum(U_total_v, 2); 
                        idx = idx + 1;

                    end
                end
            end


            U = sum(U_total, 2); % Sum across columns
            v_normal = -dot(U, NormalHere); % Normal velocity
            v_induced{iCorpo}(ChordPanel_i, SpanPanel_i) = v_normal; 
            alpha_induced{iCorpo}(ChordPanel_i, SpanPanel_i) = -atand(v_normal / U_Inf_Mag);
        end
    end
end

% Initialize variables for results
D_2D = cell(config.NCorpi, 1); % Drag distribution for each body
C_d = cell(config.NCorpi, 1); % Drag coefficient distribution for each body
C_D = zeros(config.NCorpi, 1); % Total drag coefficient for each body

% Loop through each body
for iCorpo = 1:config.NCorpi
    % Initialize drag distribution for the current body
    D_2D{iCorpo} = zeros(2 * config.SemiSpanwiseDiscr(iCorpo), 1);

    % Loop through each spanwise segment
    for k = 1:2 * config.SemiSpanwiseDiscr(iCorpo)
        % Sum contributions across chordwise panels
        for i = 1:config.ChordwiseDiscr(iCorpo)
            D_2D{iCorpo}(k) = D_2D{iCorpo}(k) + ...
                rho * U_Inf_Mag * Gamma{iCorpo}(i, k) * ...
                cos(deg2rad(config.DihedralAngle(iCorpo))) * ...
                sind(-alpha_induced{iCorpo}(i, k));
        end

        % Compute the drag coefficient for the spanwise section
        C_d{iCorpo}(k) = D_2D{iCorpo}(k) / (0.5 * rho * U_Inf_Mag^2);
    end
end

D_compl = zeros(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
    % Compute total drag for the current body
    Deltab = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo); % Panel spanwise width

    for i = 1:2 * config.SemiSpanwiseDiscr(iCorpo)
        D_compl(iCorpo) = D_compl(iCorpo) + Deltab * D_2D{iCorpo}(i);
    end

    % Compute total drag coefficient for the current body
    C_D(iCorpo) = D_compl(iCorpo) / (0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));
end
    C_D_totale_per_angle(angleIndex) = sum(C_D);
end
% CL plot
figure(1);
hold on;
plot(RotationY,C_L_totale_per_angle,'Color',[1 0.647 0],'LineWidth',2);
xlabel('\alpha (°)', 'FontSize',20,'FontWeight','bold');
ylabel('Coefficiente di Portanza','FontSize',20,'FontWeight','bold');
grid on;
set(gca,'FontSize',20,'FontWeight','normal');
legend('Cessna 172', 'Diamond DA40', 'Location','Best','FontSize',20);
hold off;

%CD plot
figure(2);
hold on;
plot(RotationY,C_D_totale_per_angle,'Color',[0.5 0 0.5],'LineWidth',2);
xlabel('\alpha (°)','FontSize',20,'FontWeight','bold');
ylabel('Coefficiente di Resistenza', 'FontSize',20,'FontWeight','bold');
grid on;
set(gca,'FontSize',20,'FontWeight','normal');
legend('Cessna 172', 'Diamond DA40', 'Location','Best', 'FontSize', 20);
hold off

% CD vs CL plot
figure(3);
hold on;
plot(C_D_totale_per_angle,C_L_totale_per_angle,'Color',[0.9,0.7,0],'LineWidth',2);
xlabel(' Coefficiente di resistenza', 'FontSize',20,'FontWeight','bold');
ylabel('Coefficiente di portanza','FontSize',20,'FontWeight','bold');
grid on;
set(gca,'FontSize',20,'FontWeight','bold');
legend('Cessna 172', 'Diamond DA40', 'Location','Best','FontSize', 20);
hold off


set(gcf, 'Color', 'w'); % imposta sfondo bianco anziché grigio
ax = gca;
ax.Clipping = 'off';
ax.XRuler.Exponent = 0; % per avere p.e. 5000 anziché 5e3 / modificate a piacimento
ax.YRuler.Exponent = 0;
% i seguenti valori credo possano essere ideali nel nostro caso
ax.LineWidth = 0.8;
ax.XAxis.FontSize = 12;  
ax.YAxis.FontSize = 12; 
exportgraphics(ax,"nome.png","Resolution",300) %  NON cambiare 300
%% Parameters for the two aircrafts

% Cessna wing
%config.RootChord = [1.625];
%config.DihedralAngle = [1.44]; 
%config.SweepAngle = [0]; 
%config.TaperRatio = [0.672]; 
%config.Span = [11];
%config.AspectRatio = [7.32]; 
%config.LEPosition_X = [0];
%config.LEPosition_Y = [0];
%config.LEPosition_Z = [0];
%config.RotationAngle_X = [0];
%config.RotationAngle_Z = [0];
%config.SemiSpanwiseDiscr = [5];
%config.ChordwiseDiscr = [5];

% Cessna tail
%config.RootChord = [1.4];
%config.DihedralAngle = [-1.44]; 
%config.SweepAngle = [0]; 
%config.TaperRatio = [0.5715]; 
%config.Span = [3.4];
%config.AspectRatio = [config.Span(2)^2/3.74]; 
%config.LEPosition_X = [7];
%config.LEPosition_Y = [0];
%config.LEPosition_Z = [-0.5];
%config.RotationAngle_X = [0];
%config.RotationAngle_Z = [0];
%config.SemiSpanwiseDiscr = [5];
%config.ChordwiseDiscr = [5];

% Diamond DA40 wing
%config.RootChord = [1.53];
%config.DihedralAngle = [5];
%config.SweepAngle = [1]; 
%config.TaperRatio = [0.647]; 
%config.Span = [11.94];
%config.AspectRatio = [10.53]; 
%config.LEPosition_X = [0];
%config.LEPosition_Y = [0];
%config.LEPosition_Z = [0];
%config.RotationAngle_X = [0];
%config.RotationAngle_Z = [0];
% Discretization options
%config.SemiSpanwiseDiscr = [5];
%config.ChordwiseDiscr = [5];

% Diamond DA40 tail
%config.RootChord = [0.7];
%config.DihedralAngle = [-5];
%config.SweepAngle = [1];
%config.TaperRatio = [0.5715]; 
%config.Span = [3.7];
%config.AspectRatio = [config.Span(2)^2/1.56]; 
%config.LEPosition_X = [6];
%config.LEPosition_Y = [0];
%config.LEPosition_Z = [-0.6];
%config.RotationAngle_X = [0];
%config.RotationAngle_Z = [0];
% Discretization options
%config.SemiSpanwiseDiscr = [5];
%config.ChordwiseDiscr = [5];


