function Color = colorPicker(colorName)
% colorPicker  Return the RGB triplet for a named plot colour.
%   Color = colorPicker(colorName) maps a colour name (e.g. 'blue', 'red3',
%   'greendark') to its RGB triplet; unknown names return black.
%
%   Inputs:
%     colorName - colour name string.
%   Outputs:
%     Color     - 1x3 RGB triplet in [0 1].

switch lower(colorName)
    % gem12
    case 'blue'
        Color = [     0,    0.4470,    0.7410];
        
    case 'blue1'
        Color = [0.435, 0.749, 0.988];
    case 'blue2'
        Color = [0.11, 0.569, 0.922];
    case 'blue3'
        Color = [0.067, 0.435, 0.722];
    case 'blue4'
        Color = [0.031, 0.278, 0.471];

    case 'red'
        Color = [0.8500,    0.3250,    0.0980];

    case 'red1'
        Color = [0.969, 0.537, 0.349];
    case 'red2'
        Color = [0.961, 0.376, 0.125];
    case 'red3'
        Color = [0.729, 0.282, 0.086];
    case 'red4'
        Color = [0.49, 0.184, 0.051];

    case 'purple1'
        Color = [0.918, 0.843, 0.98];
    case 'purple2'
        Color = [0.796, 0.576, 0.969];
    case 'purple3'
        Color = [0.659, 0.278, 0.961];
    case 'purple4'
        Color = [0.447, 0.071, 0.749];
    case 'yellow'
        Color = [0.9290,    0.6940,    0.1250];
    case 'purpledark'
        Color = [0.4940,    0.1840,    0.5560];
    case 'green'
        Color = [0.4660,    0.6740,    0.1880];

    case 'green1'
        Color = [0.58, 0.929, 0.125];
    case 'green2'
        Color = [0.486, 0.788, 0.086];
    case 'green3'
        Color = [0.322, 0.529, 0.047];
    case 'green4'
        Color = [0.216, 0.361, 0.027];

    case 'greendark'
        Color = [0.4660,    0.6740,    0.1880].*2/3; % a shade darker 
    case 'greenbright'
        Color = ([0.4660,    0.6740,    0.1880].*2+1)/3;  % a shade brighter 
    case 'cyan' 
        Color = [0.3010,    0.7450,    0.9330];
    case 'reddark'
        Color = [0.6350,    0.0780,    0.1840];
        % Color = [0.8500,    0.3250,    0.0980].*2/3;
    case 'yellowbright'
        Color = [1.0000,    0.8390,    0.0390];
    case 'bluebright'
        Color = [0.3960,    0.5090,    0.9920];
    case 'redbright'
        Color = [1.0000,    0.2700,    0.2270];
    case 'teal'
        Color = [     0,    0.6390,    0.6390];
    case 'brown'
        Color = [0.7960,    0.5170,    0.3640];

    case 'orange'
        Color = [0.9960,    0.5640,    0.2620];
    case 'pink'
        Color = [0.9490,    0.4030,    0.7720];
    case 'purple'
        Color = [0.7520,    0.3600,    0.9840];
    case 'gray'
        Color = [   0.4,       0.4,       0.4];
    case 'black'
        Color = [     0,         0,         0];
    otherwise
        Color = [     0,         0,         0];
end