

yesBL = zeros(0,2);
yesSound = zeros(0,2);
noBL = zeros(0,2);
noSound = zeros(0,2);

for u = 1:127
    if unitData(u).type~=0
        if ismember(unitData(u).neuronNumber,responsiveUnits.laserResponsiveUnits)
            bl = [unitData(u).meanResponse(2,2) unitData(u).meanResponse(2,4)];
            yesBL = [yesBL; bl];
            
            s = [unitData(u).meanResponse(1,4) unitData(u).meanResponse(3,4)];
            yesSound = [yesSound; s];
        else
            bl = [unitData(u).meanResponse(2,2) unitData(u).meanResponse(2,4)];
            noBL = [noBL; bl];
            
            s = [unitData(u).meanResponse(1,4) unitData(u).meanResponse(3,4)];
            noSound = [noSound; s];
        end
    end
end