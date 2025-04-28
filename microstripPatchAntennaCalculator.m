function microstripPatchAntennaCalculator()
    % Enhanced Microstrip Patch Antenna Calculator
    % Calculates antenna parameters using exact equations for accuracy
    
    clc;
    disp('Enhanced Microstrip Patch Antenna Calculator');
    disp('------------------------------------------');
    disp('1. Calculate dimensions and parameters from frequency');
    disp('2. Calculate frequency and parameters from dimensions');
    userChoice = input('Select option (1-2): ');
    
    % Get common parameters needed for all calculations
    substratePermittivity = input('Enter substrate dielectric constant (εr): ');
    substrateThickness = input('Enter substrate thickness (mm): ');
    
    switch userChoice
        case 1
            % Mode 1: Frequency -> Dimensions
            resonantFrequency = input('Enter resonant frequency (GHz): ');
            
            % Calculate physical dimensions
            [patchWidth, patchLength, effectivePermittivity] = ...
                calculateDimensionsFromFrequency(resonantFrequency, substratePermittivity, substrateThickness);
            
            % Calculate all other parameters
            calculatedFrequency = resonantFrequency; % This is our input
            bandwidth = calculateBandwidth(resonantFrequency, substratePermittivity, substrateThickness, patchWidth, patchLength);
            [higherModeFrequencies, extendedLength] = ...
                calculateHigherModes(resonantFrequency, effectivePermittivity, substrateThickness, patchWidth, patchLength);
            [edgeResistancePlus, edgeResistanceMinus] = ...
                calculateEdgeResistance(resonantFrequency, effectivePermittivity, patchWidth, patchLength, substrateThickness);
            
        case 2
            % Mode 2: Dimensions -> Frequency
            patchWidth = input('Enter patch width (mm): ');
            patchLength = input('Enter patch length (mm): ');
            
            % Calculate effective dielectric constant
            effectivePermittivity = calculateEffectivePermittivity(substratePermittivity, substrateThickness, patchWidth);
            
            % Calculate frequency from dimensions
            [calculatedFrequency, extendedLength] = ...
                calculateFrequencyFromDimensions(effectivePermittivity, substrateThickness, patchWidth, patchLength);
            
            % Calculate all other parameters
            bandwidth = calculateBandwidth(calculatedFrequency, substratePermittivity, substrateThickness, patchWidth, patchLength);
            [higherModeFrequencies, ~] = ...
                calculateHigherModes(calculatedFrequency, effectivePermittivity, substrateThickness, patchWidth, patchLength);
            [edgeResistancePlus, edgeResistanceMinus] = ...
                calculateEdgeResistance(calculatedFrequency, effectivePermittivity, patchWidth, patchLength, substrateThickness);
            
        otherwise
            disp('Invalid choice');
            return;
    end
    
    % Display all results
    disp(' ');
    disp('MICROSTRIP PATCH ANTENNA PARAMETERS:');
    disp('-----------------------------------');
    disp(['Resonant Frequency (TM010 mode): ' num2str(calculatedFrequency) ' GHz']);
    disp(['Effective Dielectric Constant: ' num2str(effectivePermittivity)]);
    disp(['Patch Width (W): ' num2str(patchWidth) ' mm']);
    disp(['Patch Length (L): ' num2str(patchLength) ' mm']);
    disp(['Extended Length (L*): ' num2str(extendedLength) ' mm']);
    disp(['Bandwidth: ' num2str(bandwidth) '%']);
    disp(['Edge Resistance (Rin0+): ' num2str(edgeResistancePlus) ' ohms']);
    disp(['Edge Resistance (Rin0-): ' num2str(edgeResistanceMinus) ' ohms']);
    
    % Display higher order modes
    disp(' ');
    disp('Higher Order Mode Frequencies:');
    disp(['TM020: ' num2str(higherModeFrequencies(2,1,1)) ' GHz']);
    disp(['TM110: ' num2str(higherModeFrequencies(1,1,1)) ' GHz']);
    disp(['TM011: ' num2str(higherModeFrequencies(1,1,2)) ' GHz']);
end

function [width, length, effectivePermittivity] = calculateDimensionsFromFrequency(resonantFreq, permittivity, thickness)
    % Calculate patch dimensions from resonant frequency
    % Inputs:
    %   resonantFreq - Resonant frequency in GHz
    %   permittivity - Substrate dielectric constant
    %   thickness - Substrate thickness in mm
    % Outputs:
    %   width - Patch width in mm
    %   length - Patch length in mm
    %   effectivePermittivity - Effective dielectric constant
    
    % Convert units
    freqHz = resonantFreq * 1e9; % GHz to Hz
    thicknessM = thickness * 1e-3; % mm to m
    speedOfLight = 3e8; % m/s
    
    % Calculate width using standard formula
    width = speedOfLight/(2*freqHz*sqrt((permittivity + 1)/2));
    
    % Calculate effective dielectric constant
    effectivePermittivity = calculateEffectivePermittivity(permittivity, thickness, width*1e3);
    
    % Calculate length using fringing field correction
    firstTerm = speedOfLight/(2*freqHz*sqrt(effectivePermittivity));
    numerator = (effectivePermittivity + 0.3)*(width/thicknessM + 0.264);
    denominator = (effectivePermittivity - 0.258)*(width/thicknessM + 0.8);
    secondTerm = 0.824*thicknessM*(numerator/denominator);
    length = firstTerm - secondTerm;
    
    % Convert to mm for output
    width = width * 1e3;
    length = length * 1e3;
end

function [resonantFreq, extendedLength] = calculateFrequencyFromDimensions(effectivePermittivity, thickness, width, length)
    % Calculate resonant frequency from patch dimensions
    % Inputs:
    %   effectivePermittivity - Effective dielectric constant
    %   thickness - Substrate thickness in mm
    %   width - Patch width in mm
    %   length - Patch length in mm
    % Outputs:
    %   resonantFreq - Resonant frequency in GHz
    %   extendedLength - Effective length with fringing in mm
    
    % Convert units
    thicknessM = thickness * 1e-3; % mm to m
    widthM = width * 1e-3; % mm to m
    lengthM = length * 1e-3; % mm to m
    speedOfLight = 3e8; % m/s
    
    % Calculate extended length due to fringing fields
    numerator = (effectivePermittivity + 0.3)*(widthM/thicknessM + 0.264);
    denominator = (effectivePermittivity - 0.258)*(widthM/thicknessM + 0.8);
    extendedLength = lengthM + 0.824*thicknessM*(numerator/denominator);
    
    % Calculate resonant frequency
    resonantFreq = speedOfLight/(2*extendedLength*sqrt(effectivePermittivity));
    resonantFreq = resonantFreq / 1e9; % Convert to GHz
    
    % Convert extendedLength back to mm
    extendedLength = extendedLength * 1e3;
end

function effectivePermittivity = calculateEffectivePermittivity(permittivity, thickness, width)
    % Calculate effective dielectric constant considering fringing fields
    % Inputs:
    %   permittivity - Substrate dielectric constant
    %   thickness - Substrate thickness in mm
    %   width - Patch width in mm
    % Output:
    %   effectivePermittivity - Effective dielectric constant
    
    effectivePermittivity = (permittivity + 1)/2 + ...
        (permittivity - 1)/2 * (1/sqrt(1 + 12*(thickness/width)));
end

function bandwidth = calculateBandwidth(resonantFreq, permittivity, thickness, width, length)
    % Calculate antenna bandwidth percentage
    % Inputs:
    %   resonantFreq - Resonant frequency in GHz
    %   permittivity - Substrate dielectric constant
    %   thickness - Substrate thickness in mm
    %   width - Patch width in mm
    %   length - Patch length in mm
    % Output:
    %   bandwidth - Bandwidth percentage
    
    % Convert units
    freqHz = resonantFreq * 1e9; % GHz to Hz
    thicknessM = thickness * 1e-3; % mm to m
    widthM = width * 1e-3; % mm to m
    lengthM = length * 1e-3; % mm to m
    speedOfLight = 3e8; % m/s
    wavelength = speedOfLight/freqHz;
    
    % Calculate bandwidth using standard formula
    bandwidth = 3.77 * freqHz * ((permittivity - 1)/permittivity^2) * ...
        (widthM*thicknessM/(lengthM*wavelength))^2;
    bandwidth = bandwidth / 1e6; % Convert to percentage
end

function [modeFrequencies, extendedLength] = calculateHigherModes(resonantFreq, effectivePermittivity, thickness, width, length)
    % Calculate higher order mode frequencies
    % Inputs:
    %   resonantFreq - Fundamental resonant frequency in GHz
    %   effectivePermittivity - Effective dielectric constant
    %   thickness - Substrate thickness in mm
    %   width - Patch width in mm
    %   length - Patch length in mm
    % Outputs:
    %   modeFrequencies - 3D array of mode frequencies (i,j,k) in GHz
    %   extendedLength - Effective length with fringing in mm
    
    % Convert units
    freqHz = resonantFreq * 1e9; % GHz to Hz
    thicknessM = thickness * 1e-3; % mm to m
    widthM = width * 1e-3; % mm to m
    lengthM = length * 1e-3; % mm to m
    speedOfLight = 3e8; % m/s
    
    % Calculate extended length due to fringing fields
    numerator = (effectivePermittivity + 0.3)*(widthM/thicknessM + 0.264);
    denominator = (effectivePermittivity - 0.258)*(widthM/thicknessM + 0.8);
    extendedLength = lengthM + 0.824*thicknessM*(numerator/denominator);
    
    % Initialize mode frequency matrix (i,j,k indices)
    modeFrequencies = zeros(2,2,2);
    
    % Calculate higher order modes (i,j,k)
    for i = 1:2
        for j = 1:2
            for k = 1:2
                modeFrequencies(i,j,k) = speedOfLight/(2*sqrt(effectivePermittivity)) * ...
                    sqrt((i/extendedLength)^2 + (j/widthM)^2 + (k/thicknessM)^2);
                modeFrequencies(i,j,k) = modeFrequencies(i,j,k) / 1e9; % Convert to GHz
            end
        end
    end
    
    % Convert extendedLength back to mm
    extendedLength = extendedLength * 1e3;
end

function [resistancePlus, resistanceMinus] = calculateEdgeResistance(resonantFreq, effectivePermittivity, width, length, thickness)
    % Calculate edge resistance of microstrip patch antenna
    % Inputs:
    %   resonantFreq - Resonant frequency in GHz
    %   effectivePermittivity - Effective dielectric constant
    %   width - Patch width in mm
    %   length - Patch length in mm
    %   thickness - Substrate height in mm
    % Outputs:
    %   resistancePlus - Edge resistance for + case (ohms)
    %   resistanceMinus - Edge resistance for - case (ohms)
    
    % Convert units
    freqHz = resonantFreq * 1e9; % GHz to Hz
    widthM = width * 1e-3; % mm to m
    lengthM = length * 1e-3; % mm to m
    thicknessM = thickness * 1e-3; % mm to m
    speedOfLight = 3e8; % m/s
    
    % Calculate wave parameters
    wavelength = speedOfLight/freqHz;
    waveNumber = 2*pi/wavelength;
    
    % Calculate conductance terms
    [selfConductance, mutualConductance] = ...
        calculateConductanceIntegrals(widthM, lengthM, waveNumber);
    
    % Calculate edge resistances
    resistancePlus = 1/(2*(selfConductance + mutualConductance));
    resistanceMinus = 1/(2*(selfConductance - mutualConductance));
end

function [selfConductance, mutualConductance] = calculateConductanceIntegrals(width, length, waveNumber)
    % Numerical integration for conductance calculations
    % Inputs:
    %   width - Patch width in meters
    %   length - Patch length in meters
    %   waveNumber - Wave number in free space
    % Outputs:
    %   selfConductance - Self conductance term
    %   mutualConductance - Mutual conductance term
    
    % Create angular sampling points (0 to 180 degrees)
    thetaDegrees = 0:1:180; 
    thetaRadians = thetaDegrees.*pi/180;
    
    % Calculate argument for sinc function
    sincArgument = cos(thetaRadians).*(waveNumber*width/2);
    
    % Calculate integrand for self conductance
    selfIntegrand = sinc(sincArgument./pi).^2 .* sin(thetaRadians).^2 .* ...
        sin(thetaRadians) .* ((pi/180)*(waveNumber*width/2)^2);
    
    % Calculate integrand for mutual conductance
    mutualIntegrand = sinc(sincArgument./pi).^2 .* sin(thetaRadians).^2 .* ...
        besselj(0, sin(thetaRadians).*(waveNumber*length)) .* ...
        sin(thetaRadians) .* ((pi/180)*(waveNumber*width/2)^2);
    
    % Sum to approximate integral
    selfConductance = sum(selfIntegrand)./(120*pi^2); 
    mutualConductance = sum(mutualIntegrand)./(120*pi^2);
end
