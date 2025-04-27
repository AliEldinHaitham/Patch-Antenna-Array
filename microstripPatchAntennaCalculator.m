function microstripPatchAntennaCalculator()
    % Enhanced Microstrip Patch Antenna Calculator
    % Using the exact equations provided for more accurate calculations
    
    clc;
    disp('Enhanced Microstrip Patch Antenna Calculator');
    disp('------------------------------------------');
    disp('1. Calculate dimensions and parameters from frequency');
    disp('2. Calculate frequency and parameters from dimensions');
    choice = input('Select option (1-2): ');
    
    % Common parameters needed for all calculations
    epsilon_r = input('Enter substrate dielectric constant (εr): ');
    h = input('Enter substrate thickness (mm): ');
    
    switch choice
        case 1
            % Mode 1: Frequency -> Dimensions
            f0 = input('Enter resonant frequency (GHz): ');
            
            % Calculate dimensions
            [W, L, epsilon_eff] = calculateDimensionsFromFreq(f0, epsilon_r, h);
            
            % Calculate all other parameters
            calculated_freq = f0; % This is our input
            BW = calculateBandwidth(f0, epsilon_r, h, W, L);
            [f_ijk, L_star] = calculateHigherModes(f0, epsilon_eff, h, W, L);
            [Rin0P, Rin0M] = calculateEdgeResistance(f0, epsilon_eff, W, L, h);
            
        case 2
            % Mode 2: Dimensions -> Frequency
            W = input('Enter patch width (mm): ');
            L = input('Enter patch length (mm): ');
            
            % Calculate effective dielectric constant
            epsilon_eff = calculateEpsilonEff(epsilon_r, h, W);
            
            % Calculate frequency from dimensions
            [calculated_freq, L_star] = calculateFrequencyFromDims(epsilon_eff, h, W, L);
            
            % Calculate all other parameters
            BW = calculateBandwidth(calculated_freq, epsilon_r, h, W, L);
            [f_ijk, ~] = calculateHigherModes(calculated_freq, epsilon_eff, h, W, L);
            [Rin0P, Rin0M] = calculateEdgeResistance(calculated_freq, epsilon_eff, W, L, h);
            
        otherwise
            disp('Invalid choice');
            return;
    end
    
    % Display all results
    disp(' ');
    disp('MICROSTRIP PATCH ANTENNA PARAMETERS:');
    disp('-----------------------------------');
    disp(['Resonant Frequency (TM010 mode): ' num2str(calculated_freq) ' GHz']);
    disp(['Effective Dielectric Constant: ' num2str(epsilon_eff)]);
    disp(['Patch Width (W): ' num2str(W) ' mm']);
    disp(['Patch Length (L): ' num2str(L) ' mm']);
    disp(['Extended Length (L*): ' num2str(L_star) ' mm']);
    disp(['Bandwidth: ' num2str(BW) '%']);
    disp(['Edge Resistance (Rin0+): ' num2str(Rin0P) ' ohms']);
    disp(['Edge Resistance (Rin0-): ' num2str(Rin0M) ' ohms']);
    
    % Display higher order modes
    disp(' ');
    disp('Higher Order Mode Frequencies:');
    disp(['TM020: ' num2str(f_ijk(2,1,1)) ' GHz']);
    disp(['TM110: ' num2str(f_ijk(1,1,1)) ' GHz']);
    disp(['TM011: ' num2str(f_ijk(1,1,2)) ' GHz']);
end

function [W, L, epsilon_eff] = calculateDimensionsFromFreq(f0, epsilon_r, h)
    % Convert units
    f0 = f0 * 1e9; % GHz to Hz
    h = h * 1e-3;  % mm to m
    
    c = 3e8; % speed of light
    
    % Calculate width using provided equation
    W = c/(2*f0*sqrt((epsilon_r + 1)/2));
    
    % Calculate effective dielectric constant
    epsilon_eff = calculateEpsilonEff(epsilon_r, h*1e3, W*1e3); % Convert to mm for function
    
    % Calculate length using provided equation
    term1 = c/(2*f0*sqrt(epsilon_eff));
    term2_num = (epsilon_eff + 0.3)*(W/h + 0.264);
    term2_den = (epsilon_eff - 0.258)*(W/h + 0.8);
    term2 = 0.824*h*(term2_num/term2_den);
    L = term1 - term2;
    
    % Convert to mm for output
    W = W * 1e3;
    L = L * 1e3;
end

function [f0, L_star] = calculateFrequencyFromDims(epsilon_eff, h, W, L)
    % Convert units
    h = h * 1e-3;  % mm to m
    W = W * 1e-3;   % mm to m
    L = L * 1e-3;   % mm to m
    
    c = 3e8; % speed of light
    
    % Calculate extended length L*
    term_num = (epsilon_eff + 0.3)*(W/h + 0.264);
    term_den = (epsilon_eff - 0.258)*(W/h + 0.8);
    L_star = L + 0.824*h*(term_num/term_den);
    
    % Calculate resonant frequency
    f0 = c/(2*L_star*sqrt(epsilon_eff));
    f0 = f0 / 1e9; % Convert to GHz
    
    % Convert L_star back to mm
    L_star = L_star * 1e3;
end

function epsilon_eff = calculateEpsilonEff(epsilon_r, h, W)
    % Calculate effective dielectric constant using provided equation
    epsilon_eff = (epsilon_r + 1)/2 + (epsilon_r - 1)/2 * (1/sqrt(1 + 12*(h/W)));
end

function BW = calculateBandwidth(f0, epsilon_r, h, W, L)
    % Convert units
    f0 = f0 * 1e9; % GHz to Hz
    h = h * 1e-3;  % mm to m
    W = W * 1e-3;  % mm to m
    L = L * 1e-3;  % mm to m
    c = 3e8;
    lambda = c/f0;
    
    % Calculate bandwidth using provided equation
    BW = 3.77 * f0 * ((epsilon_r - 1)/epsilon_r^2) * (W*h/(L*lambda))^2;
    BW = BW / 1e6; % Convert to percentage
end

function [f_ijk, L_star] = calculateHigherModes(f0, epsilon_eff, h, W, L)
    % Convert units
    f0 = f0 * 1e9; % GHz to Hz
    h = h * 1e-3;  % mm to m
    W = W * 1e-3;  % mm to m
    L = L * 1e-3;  % mm to m
    c = 3e8;
    
    % Calculate extended length L*
    term_num = (epsilon_eff + 0.3)*(W/h + 0.264);
    term_den = (epsilon_eff - 0.258)*(W/h + 0.8);
    L_star = L + 0.824*h*(term_num/term_den);
    
    % Initialize mode frequency matrix
    f_ijk = zeros(2,2,2);
    
    % Calculate higher order modes (i,j,k)
    for i = 1:2
        for j = 1:2
            for k = 1:2
                f_ijk(i,j,k) = c/(2*sqrt(epsilon_eff)) * sqrt((i/L_star)^2 + (j/W)^2 + (k/h)^2);
                f_ijk(i,j,k) = f_ijk(i,j,k) / 1e9; % Convert to GHz
            end
        end
    end
    
    % Convert L_star back to mm
    L_star = L_star * 1e3;
end

function [Rin0P, Rin0M] = calculateEdgeResistance(f0, epsilon_eff, W, L, h)
    % Calculate edge resistance of microstrip patch antenna
    % Inputs:
    %   f0 - Resonant frequency (GHz)
    %   epsilon_eff - Effective dielectric constant
    %   W - Patch width (mm)
    %   L - Patch length (mm)
    %   h - Substrate height (mm)
    % Outputs:
    %   Rin0P - Edge resistance for + case (ohms)
    %   Rin0M - Edge resistance for - case (ohms)
    
    % Convert units
    f0 = f0 * 1e9;  % GHz to Hz
    W = W * 1e-3;   % mm to m
    L = L * 1e-3;   % mm to m
    h = h * 1e-3;   % mm to m
    
    c = 3e8; % speed of light
    lambda_o = c/f0;
    ko = 2*pi/lambda_o;
    
    % Calculate conductance terms
    [G1, G12] = sintegr(W, L, ko);
    
    % Calculate edge resistances
    Rin0P = 1/(2*(G1 + G12));
    Rin0M = 1/(2*(G1 - G12));
end

function [G1, G12] = sintegr(W, L, ko)
    % Numerical integration for conductance calculations
    % Inputs:
    %   W - Patch width (m)
    %   L - Patch length (m)
    %   ko - Wave number in free space
    % Outputs:
    %   G1 - Self conductance
    %   G12 - Mutual conductance
    
    th = 0:1:180; 
    t = th.*pi/180;
    ARG = cos(t).*(ko*W/2);
    res1 = sum(sinc(ARG./pi).^2.*sin(t).^2.*sin(t).*((pi/180)*(ko*W/2)^2));
    res12 = sum(sinc(ARG./pi).^2.*sin(t).^2.*besselj(0,sin(t).*(ko*L)).*sin(t).*((pi/180)*(ko*W/2)^2));
    G1 = res1./(120*pi^2); 
    G12 = res12./(120*pi^2);
end