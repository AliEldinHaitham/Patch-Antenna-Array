function Taylor_array()
clc;
clear all;

% USER INPUT PARAMETERS
Nelem = 0;
while Nelem < 2
    Nelem = floor(input('NUMBER OF ELEMENTS = '));
end

d = input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) = ');
SLL_dB =0;
n_bar=0;


while SLL_dB >= 0
    SLL_dB = -input('SIDE LOBE LEVEL (IN dB, e.g., 20) = ');
end

while n_bar < 1 || mod(n_bar,1) ~= 0
    n_bar = floor(input('TAYLOR n_bar: NO. OF NEARLY EQUAL SIDELOBES (>=1) = '));
end

beta = 0; % Broadside by default

theta = linspace(0, pi, 1800);

% TAYLOR WEIGHTS
I_taylor = taylorwin(Nelem, n_bar, SLL_dB); % Now includes n_bar
I_taylor = I_taylor / max(I_taylor);   % Normalize

% Compute AF over theta
AF = zeros(size(theta));
for n = 1:Nelem
    phase = 2 * pi * d * (n - (Nelem + 1)/2) * cos(theta) + beta;
    AF = AF + I_taylor(n) * exp(1j * phase);
end
AF = abs(AF);
AF = AF / max(AF); % Normalize

U = AF.^2;
dtheta = pi / length(theta);
Prad = 2 * pi * sum(U .* sin(theta) .* dtheta);
D = 4 * pi * U ./ (Prad + eps);
DdB = 10 * log10(D + eps);

% Find HPBW
[hp, ~] = hpbw(AF);
Do = max(D);
DodB = max(DdB);

% =============================
% FIGURE 1: ARRAY FACTOR (dB)
% =============================
figure;
plot(theta*180/pi, 20*log10(AF), 'm', 'LineWidth', 2);
xlabel('\theta (degrees)');
ylabel('Normalized |AF| (dB)');
grid on;
axis([0 180 -60 0]);
title(sprintf('Taylor Array Pattern (N=%d, SLL=%.1fdB, bar(n)=%d)', Nelem, SLL_dB, n_bar), ...
    'Interpreter','latex');
text(180, -5, ['HPBW = ', num2str(max(hp)), ' deg'], ...
    'HorizontalAlignment','right');

% =============================
% ADD MARKERS AT Y = -3 dB
% =============================
AF_dB = 20 * log10(AF);
threshold = -3;

% Find indices where AF crosses -3 dB
cross_indices = find(diff(sign(AF_dB - threshold)) ~= 0);

if ~isempty(cross_indices)
    hold on;
    for k = 1:length(cross_indices)
        idx = cross_indices(k);

        % Interpolate to get more accurate crossing point
        x_before = theta(idx) * 180/pi;
        x_after = theta(idx+1) * 180/pi;
        y_before = AF_dB(idx);
        y_after = AF_dB(idx+1);

        % Linear interpolation
        slope = (y_after - y_before) / (x_after - x_before);
        if slope == 0
            continue;
        end
        x_cross = x_before + (threshold - y_before)/slope;

        % Plot vertical line at crossing point
        plot([x_cross x_cross], [min(ylim) threshold], 'k--', 'LineWidth', 1);

        % Plot a marker at the intersection point
        plot(x_cross, threshold, 'ro', 'MarkerFaceColor', 'r');
    end
    hold off;
end
% =============================
% FIGURE 2: ARRAY FACTOR (Linear)
% =============================
figure;
plot(theta*180/pi, AF, 'm', 'LineWidth', 2);
xlabel('\theta (degrees)');
ylabel('Normalized |AF|');
grid on;
axis([0 180 0 1]);

% =============================
% DIRECTIVITY PLOTS
% =============================
figure;
subplot(2,1,1);
plot(theta*180/pi, D, 'r', 'LineWidth', 2);
ylabel('Directivity (dimensionless)');
grid on;
title('Directivity vs \theta', 'FontSize', 12);
text(180, 1.5*Do, ['D_0 = ', num2str(Do)], 'HorizontalAlignment','right');

subplot(2,1,2);
plot(theta*180/pi, DdB, 'b', 'LineWidth', 2);
xlabel('\theta (degrees)');
ylabel('Directivity (dB)');
grid on;
text(180, 1.5*DodB, ['D_0 = ', num2str(DodB), ' dB'], 'HorizontalAlignment','right');
ylim([SLL_dB - 30, DodB]); % Ensure it's [min max]

% =============================
% FIGURE 5: EXCITATION COEFFICIENTS
% =============================
x = (1:Nelem) - (Nelem+1)/2; % Element positions
figure;
[AX, H1, H2] = plotyy(x, I_taylor, x, angle(I_taylor)*180/pi);
set(get(AX(1),'Ylabel'),'String','Amplitude','Color','r');
set(get(AX(2),'Ylabel'),'String','Phase (deg)','Color','b');
set(H1, 'LineStyle','-','Marker','o','Color','r','LineWidth',1.5);
set(H2, 'LineStyle',':','Marker','s','Color','b','LineWidth',1.5);
legend('Amplitude','Phase');
xlabel('Element Index');
title('Excitation Amplitude and Phase', 'FontSize', 12);

end % End main function

% ---------------------
% HPBW Calculation
% ---------------------
function[hp,thmax] = hpbw(AF)
tol = 0.001;
imax = 0;
j = 0;
M = length(AF);
root = [];
hp = [];

for i = 2:M
    if AF(i) > AF(i-1) && AF(i) > (1 - tol)
        imax = imax + 1;
        thmax(imax) = i;
    end
    if i > 1
        y1 = AF(i) - 0.707;
        y2 = AF(i-1) - 0.707;
        if y1 * y2 < 0
            j = j + 1;
            root(j) = (i - 1) + (-y2) * (1)/(y1 - y2);
            if j >= 2
                hp(end+1) = root(j) - root(j-1);
            end
        end
    end
end
if isempty(hp)
    hp = 0;
else
    hp = hp * 180 / M; % Convert index to degrees
end
end
