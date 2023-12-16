clear all; close all; clc;

% Constants
P_LED = 2; % transmitted optical power by individual LED
Adet = 1 * 10 ^ -4; % detector physical area of a PD
Ts = 1; % gain of an optical filter
index = 1.5; % refractive index of a lens at a PD
FOV = 140; % FOV of a receiver
G_Con = (index ^ 2) / (sind(FOV) .^ 2); % gain of an optical concentrator
lx = 2; ly = 2; h = 3; % room dimensions in meters
XT = 0; YT = 0; % Position of the single LED
resolution = 20;
Nx = lx * resolution; Ny = ly * resolution;
x = linspace(-lx / 2, lx / 2, Nx);
y = linspace(-ly / 2, ly / 2, Ny);
[XR, YR] = meshgrid(x, y);

% Prepare for animation
filenameWatt = '3D_power_curve_Watt.gif';
filenameDBm = '3D_power_curve_dBm.gif';

power_at_corner = [];
power_at_center = [];

for ml = 1:20 % Lambertian order range

    D = sqrt((XR - XT) .^ 2 + (YR - YT) .^ 2 + h ^ 2); % Distance vector
    cosphi = h ./ D; % Angle vector
    receiver_angle = acosd(cosphi);
    H = (ml + 1) * Adet .* cosphi .^ (ml + 1) ./ (2 * pi .* D .^ 2); % Channel DC gain
    P_rec = P_LED .* H .* Ts .* G_Con; % Received power
    P_rec(abs(receiver_angle) > FOV) = 0;
    P_rec_dBm = 10 * log10(P_rec * 1000);

    % Find the received power at the center and one of the corners
    P_center = P_rec(Ny/2, Nx/2); % Center of the room
    P_corner = P_rec(1, 1); % Top-left corner of the room

    % Convert powers to dBm for comparison
    P_center_dBm = 10 * log10(P_center * 1000);
    P_corner_dBm = 10 * log10(P_corner * 1000);

    % Store the powers for each Lambertian order
    power_at_center = [power_at_center P_center_dBm];
    power_at_corner = [power_at_corner P_corner_dBm];

    % Plot in Watt
    figure(1)
    surfc(x, y, P_rec);
    xlabel('Room Length (m)');
    ylabel('Room Width (m)');
    zlabel('Power (Watt)');
    title(['Received Power (Watt) - Lambertian = ', num2str(ml)]);
    colorbar
    drawnow
    % Save as PNG and SVG
    saveas(gcf, ['3D_power_curve_Watt_', num2str(ml), '.png']);
    saveas(gcf, ['3D_power_curve_Watt_', num2str(ml), '.svg']);
    % Append to GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if ml == 1
        imwrite(imind, cm, filenameWatt, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, filenameWatt, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end

    % Plot in dBm
    figure(2)
    surfc(x, y, P_rec_dBm);
    xlabel('Room Length (m)');
    ylabel('Room Width (m)');
    zlabel('Power (dBm)');
    title(['Received Power (dBm) - Lambertian = ', num2str(ml)]);
    colorbar
    drawnow
    % Save as PNG and SVG
    saveas(gcf, ['3D_power_curve_dBm_', num2str(ml), '.png']);
    saveas(gcf, ['3D_power_curve_dBm_', num2str(ml), '.svg']);
    % Append to GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if ml == 1
        imwrite(imind, cm, filenameDBm, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, filenameDBm, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end

% Find the difference that close to 6 dB between the center and the corner
power_difference = power_at_center - power_at_corner;
% Find the Lambertian order that gives the closest to 6 dB difference
[~, index] = min(abs(power_difference - 6));

fprintf('The Lambertian order that gives the closest to 6 dB difference is %d\n', index);
fprintf('The received power at the center is %f dBm\n', power_at_center(index));
fprintf('The received power at the corner is %f dBm\n', power_at_corner(index));
fprintf('The difference between the center and the corner is %f dB\n', power_difference(index));