fileID = fopen('MUOS.42','r');
formatSpec = '%f';
sizeMUOS = [6 Inf];
MUOS = fscanf(fileID, formatSpec, sizeMUOS);

fileID = fopen('Orbcomm.42','r');
formatSpec = '%f';
sizeOrbcomm = [6 Inf];
Orbcomm = fscanf(fileID, formatSpec, sizeOrbcomm);

StartEpoch = MUOS(1,1);
DispDuration = 60*60*24* 8; % 8 days
EndEpoch = StartEpoch + DispDuration;

[row, col] = find(MUOS==EndEpoch);
LastColMUOS = col(length(col));

SampleMUOS = MUOS(1:sizeMUOS(1), 1:LastColMUOS);

[row, col] = find(Orbcomm==EndEpoch);
LastColOrbcomm = col(length(col));

SampleOrbcomm = Orbcomm(1:sizeOrbcomm(1), 1:LastColOrbcomm);

% MUOS

lon_MUOS = SampleMUOS(2, 1:LastColMUOS);
lat_MUOS = SampleMUOS(3, 1:LastColMUOS);
SNR_D_MUOS = SampleMUOS(4, 1:LastColMUOS);
SNR_R_MUOS = SampleMUOS(5, 1:LastColMUOS);

% Orbcomm

lon_Orbcomm = SampleOrbcomm(2, 1:LastColOrbcomm);
lat_Orbcomm = SampleOrbcomm(3, 1:LastColOrbcomm);
SNR_D_Orbcomm = SampleOrbcomm(4, 1:LastColOrbcomm);
SNR_R_Orbcomm = SampleOrbcomm(5, 1:LastColOrbcomm);

% Plot

scatter3(lon_MUOS, lat_MUOS, SNR_D_MUOS);
