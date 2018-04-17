% landslide Calculate probability of a landslide given a ShakeMap and Landslide data files.
% lsprob = landslide(xmlfile,slopefiles,cohesion,friction);
% Input:
%  - xmlfile     ShakeMap grid XML file name.
%  - slopefiles  Cell array of landslide input slope angle files (in any order).
%  - cohesion    Soil/rock cohesion data file name.
%  - friction    Internal friction angle data file name.
% Output:
%  - lsprob      A grid structure where each cell of the grid contains landslide
%                probability of 1,10,30,50,70,90, or 99 percent.
function lsprob = landslide(xmlfile,slopefiles,cohesion,friction)

    GAMMA = 15.7; % soil unit weight (kN/m^3)
    ZETA = 2.4; %slope-normal failure depth (m)
    g = 9.8; %gravitational acceleration (m/s^2)
    C1 = 0.215; %additive constant in newmark displacement calculation
    C2 = 2.341; %first exponential constant
    C3 = -1.438; %second exponential constant
    NODATA_COHESION = 1;
    NODATA_FRICTION = 26;
    DNTHRESH = 5; %Newmark displacement threshold (cm)
    FSTHRESH = 1.01; %Safety factor threshold
    ACTHRESH = 0.05; %Critical acceleration threshold
    [smgrid,event] = readsmgrid(xmlfile);
    [nrows,ncols,nbands] = size(smgrid.grid);
    bounds = getgridbounds(smgrid);
    for i=1:length(slopefiles)
        gridfile = slopefiles{i};
        grids(i) = readesri(gridfile,bounds);
    end
    
    cogrid = readesri(cohesion,bounds);
    frgrid = readesri(friction,bounds);
    
    nco = length(find(~isnan(cogrid.grid(:))));
    nfr = length(find(~isnan(cogrid.grid(:))));
    if nco && nfr
        %Adjust cohesion down by a factor of 10
        cogrid.grid = cogrid.grid./10;
    else
        [m,n] = size(cogrid.grid);
        cogrid.grid = ones(m,n).*NODATA_COHESION;
        frgrid.grid = ones(m,n).*NODATA_FRICTION;
    end
    
    
   
    %resample the shakemap to the slope grids
    smgrid = samplegeogrid(smgrid,grids(1),'linear');
    ipga = find(strcmpi(event.bandnames,'pga'));
    if isempty(ipga)
        fprintf('Could not find PGA layer!  Returning.\n');
        lsprob = [];
        return;
    end
    pga = smgrid.grid(:,:,ipga)./100;
    
    binarysum = zeros(size(cogrid.grid));
    for i=1:length(grids)
        slope = grids(i);
        slope.grid = slope.grid./100;
        FS = cogrid.grid./(GAMMA.*ZETA.*sind(slope.grid)) + (tand(frgrid.grid)./tand(slope.grid));
        FS(FS <= FSTHRESH) = FSTHRESH;
        ac = (FS - 1) .* sind(slope.grid);
        ac(ac <= ACTHRESH) = ACTHRESH;
        p1 = (1 - ac./pga).^C2;
        p2 = (ac./pga).^C3;
        logdn = C1 + log(p1 .* p2);
        dn = exp(logdn);
        dn(pga < ac) = 1;
        iunstable = find(dn >= DNTHRESH);
        istable = find(dn < DNTHRESH);
        dn(iunstable) = 1;
        dn(istable) = 0;
        binarysum = binarysum + dn;
    end
    
    binarysum(binarysum == 1) = 0.01;
    binarysum(binarysum == 2) = 0.10;
    binarysum(binarysum == 3) = 0.30;
    binarysum(binarysum == 4) = 0.50;
    binarysum(binarysum == 5) = 0.70;
    binarysum(binarysum == 6) = 0.90;
    binarysum(binarysum == 7) = 0.99;
    binarysum(isnan(binarysum)) = 0;
    lsprob = struct();
    lsprob.grid = binarysum;
    lsprob.ulxmap = smgrid.ulxmap;
    lsprob.ulymap = smgrid.ulymap;
    lsprob.xdim = smgrid.xdim;
    lsprob.ydim = smgrid.ydim;
    lsprob.nbands = 1;
    lsprob.bandnames = {'Landslide probability'};