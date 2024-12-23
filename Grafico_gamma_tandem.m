close all
clear all
clc
%% Test Case 1

U_Inf_Mag = 60;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;

config.NCorpi = 2;

config.RootChord = [1.625,1.4];
config.DihedralAngle = [1.44,-1.44]; % [°]
config.SweepAngle = [0,0]; % [°]
config.TaperRatio = [0.672,0.5715]; 
config.Span = [11,3.4];
config.AspectRatio = [7.32, config.Span(2)^2 /3.74]; 
config.LEPosition_X = [0,7];
config.LEPosition_Y = [0,0];
config.LEPosition_Z = [0,-0.5];

config.RotationAngle_X = [0];
config.RotationAngle_Y = [1.5,-1.5];
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [5,5];
config.ChordwiseDiscr = [5,5];

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

nPanelsCorpo1 = 2 * config.SemiSpanwiseDiscr(1);
nPanelsCorpo2 = 2 * config.SemiSpanwiseDiscr(2);

span_corpo1 = linspace(-config.SemiSpan(1), config.SemiSpan(1), nPanelsCorpo1);
span_corpo2 = linspace(-config.SemiSpan(2), config.SemiSpan(2), nPanelsCorpo2);

Gamma_totale_corpo1 = Gamma{1};
Gamma_totale_corpo2 = Gamma{2};

Gamma_interpolata_corpo2 = zeros(config.ChordwiseDiscr(1), 2 * config.SemiSpanwiseDiscr(1));
for k = 1:config.ChordwiseDiscr(1)
    for j = 1:2 * config.SemiSpanwiseDiscr(1)
        if span_corpo1(j) >= -config.SemiSpan(2) && span_corpo1(j) <= config.SemiSpan(2)
    Gamma_interpolata_corpo2(k,j) = spline(span_corpo2, Gamma_totale_corpo2(k,:), span_corpo1(j));
        else
            Gamma_interpolata_corpo2(k,j) = 0;
        end
    end
end
Gamma_totale_combined = Gamma_totale_corpo1 + Gamma_interpolata_corpo2;
Gamma_totale_combined_somma = sum(Gamma_totale_combined,1);
%Grafico della distribuzione di circolazione lungo lo span
figure(1);
hold on
plot(span_corpo1, Gamma_totale_combined_somma, '-', 'LineWidth',2, 'Color',[0.9,0.7,0]);
xlabel('Posizione lungo lo span, y','FontSize',12,'FontWeight','bold');
ylabel('Circolazione, \Gamma','FontSize',12,'FontWeight','bold');
title('Distribuzione di Circolazione lungo lo Span','FontSize',14,'FontWeight','bold');
grid on;
set(gca,'FontSize',10,'FontWeight','normal');
hold off



U_Inf_Mag = 60;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;

config.NCorpi = 2;

config.RootChord = [1.53,0.7];
config.DihedralAngle = [5,-5]; % [°]
config.SweepAngle = [1,1]; % [°]
config.TaperRatio = [0.647,0.5715]; 
config.Span = [11.94,3.7];
config.AspectRatio = [10.53, config.Span(2)^2 /1.56]; 
config.LEPosition_X = [0,6];
config.LEPosition_Y = [0,0];
config.LEPosition_Z = [0,0.6];

config.RotationAngle_X = [0];
config.RotationAngle_Y = [1.5,-1.5];
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [5,5];
config.ChordwiseDiscr = [5,5];

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

nPanelsCorpo1 = 2 * config.SemiSpanwiseDiscr(1);
nPanelsCorpo2 = 2 * config.SemiSpanwiseDiscr(2);

span_corpo1 = linspace(-config.SemiSpan(1), config.SemiSpan(1), nPanelsCorpo1);
span_corpo2 = linspace(-config.SemiSpan(2), config.SemiSpan(2), nPanelsCorpo2);

Gamma_totale_corpo1 = Gamma{1};
Gamma_totale_corpo2 = Gamma{2};

Gamma_interpolata_corpo2 = zeros(config.ChordwiseDiscr(1), 2 * config.SemiSpanwiseDiscr(1));
for k = 1:config.ChordwiseDiscr(1)
    for j = 1:2 * config.SemiSpanwiseDiscr(1)
        if span_corpo1(j) >= -config.SemiSpan(2) && span_corpo1(j) <= config.SemiSpan(2)
    Gamma_interpolata_corpo2(k,j) = spline(span_corpo2, Gamma_totale_corpo2(k,:), span_corpo1(j));
        else
            Gamma_interpolata_corpo2(k,j) = 0;
        end
    end
end
Gamma_totale_combined = Gamma_totale_corpo1 + Gamma_interpolata_corpo2;
Gamma_totale_combined_somma = sum(Gamma_totale_combined,1);
 %Grafico della distribuzione di circolazione lungo lo span
figure(1);
hold on
plot(span_corpo1, Gamma_totale_combined_somma, '-', 'LineWidth',2, 'Color',[0.2,0.6,0.2]);
xlabel('Posizione lungo lo span, y','FontSize',30,'FontWeight','bold');
ylabel('Circolazione, \Gamma','FontSize',30,'FontWeight','bold');
title('Distribuzione di Circolazione lungo lo Span','FontSize',30,'FontWeight','bold');
grid on;
set(gca,'FontSize',10,'FontWeight','normal');
legend('Cessna 172', 'Diamond DA40', 'Location', 'Best', 'FontSize', 20)
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