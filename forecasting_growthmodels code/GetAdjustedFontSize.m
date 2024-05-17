
function adjustedFontSize=GetAdjustedFontSize()

% Retrieve screen DPI
screenSize = get(0, 'ScreenSize');
screenDPI = get(0, 'ScreenPixelsPerInch');

% Define a base font size (this is a font size that looks good on a standard 96 DPI screen)
baseFontSize = 24;

% Calculate scaling factor based on DPI (assuming base DPI is 96)
baseDPI = 96;
scalingFactor = screenDPI / baseDPI;

% Calculate adjusted font size
adjustedFontSize = min(baseFontSize * scalingFactor,24);
