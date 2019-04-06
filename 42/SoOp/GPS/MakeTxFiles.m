
input_Filename = 'TLEinfo.txt';
inputFile = fopen(input_Filename, 'r');

% Read TLE input file

while ~feof(inputFile)
%% Orbcomm
%      FMnum = fscanf(inputFile,'ORBCOMM FM%d\n',1);
%      TLE1 = fgets(inputFile,512);
%      TLE2 = fgets(inputFile,512);
% 
% %     TLE_Filename = sprintf('TLE_FM%d.txt',FMnum);
%     
%     Orb_Filename = sprintf('Orb_FM%d.txt',FMnum);
%     system(['cp Orb_FM.txt ' Orb_Filename]);
% 
%     SC_Filename = sprintf('SC_FM%d.txt',FMnum);
%     system(['cp SC_FM.txt ' SC_Filename]);
%     
% %     TLEfile = fopen(TLE_Filename, 'w');
% %     fprintf(TLEfile,'FM%d\n',FMnum);
% %     fprintf(TLEfile,'%s',TLE1);
% %     fprintf(TLEfile,'%s',TLE2);
% %     fclose(TLEfile);
%     
%     s1 = sprintf('ORBCOMM FM%d	              !  Description\n', FMnum);
%     s2 = sprintf('TLE  "ORBCOMM FM%02d"          !  TLE or TRV format, Label to find in file\n', FMnum);
%     s3 = sprintf('"Inp_TLE.txt"               !  File name\n');
%     OverwriteLineInFile(Orb_Filename,2,s1);
%     OverwriteLineInFile(Orb_Filename,23,s2);
%     OverwriteLineInFile(Orb_Filename,24,s3);
% 
%     s2 = sprintf('"FM%d"                        !  Label\n', FMnum);
%     OverwriteLineInFile(SC_Filename,2,s1);
%     OverwriteLineInFile(SC_Filename,3,s2);
%     
% %     nextline = fgets(inputFile,512);
% %     if ~ischar(nextline)
% %         break;
% %     end
%% GPS
     PRNnum = fscanf(inputFile,'PRN %d\n',1);
     TLE1 = fgets(inputFile,512);
     TLE2 = fgets(inputFile,512);
    
    Orb_Filename = sprintf('Orb_PRN%d.txt',PRNnum);
    system(['cp Orb.txt ' Orb_Filename]);

    SC_Filename = sprintf('SC_PRN%d.txt',PRNnum);
    system(['cp SC.txt ' SC_Filename]);
    
    s1 = sprintf('GPS PRN %d	              !  Description\n', PRNnum);
    s2 = sprintf('TLE  "PRN %02d"          !  TLE or TRV format, Label to find in file\n', PRNnum);
    s3 = sprintf('"Inp_TLE.txt"               !  File name\n');
    OverwriteLineInFile(Orb_Filename,2,s1);
    OverwriteLineInFile(Orb_Filename,23,s2);
    OverwriteLineInFile(Orb_Filename,24,s3);

    s2 = sprintf('"PRN%d"                        !  Label\n', PRNnum);
    s3 = sprintf('GenTx3SpriteAlpha.ppm         !  Sprite File Name\n');

    OverwriteLineInFile(SC_Filename,2,s1);
    OverwriteLineInFile(SC_Filename,3,s2);
    OverwriteLineInFile(SC_Filename,4,s3);
end

fclose(inputFile);