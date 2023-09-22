function data = extractSIRdataColumns(SIRData)
    
    data.lDosages = [-8.96580000000000	-7.96580000000000	-6.96580000000000	-5.96580000000000	-5.05890000000000	-4.05890000000000	-3	-2	-1	0	1	2	3	4	5	6	7	8	9];
    data.lDosagesString = ['-8.9658'	'-7.9658'	'-6.9658'	'-5.9658'	'-5.0589'	'-4.0589'	'-3'	'-2'	'-1'	'0'	'1'	'2'	'3'	'4'	'5'	'6'	'7'	'8'	'9'];
    
    drugs = SIRData(3:end,1);
    data.drugs = {drugs{:}};
    bugs = SIRData(3:end,2);
    data.bugs = {bugs{:}};

    data.Nclusters = cell2mat(SIRData(3:end,3));
    data.GMMconverged = cell2mat(SIRData(3:end,4));
    data.GMMRsquared = cell2mat(SIRData(3:end,5));

    data.GMMSmean = cell2mat(SIRData(3:end,6));
    data.GMMSsigma = cell2mat(SIRData(3:end,7));
    data.GMMImean = cell2mat(SIRData(3:end,8));
    data.GMMIsigma = cell2mat(SIRData(3:end,9));
    data.GMMRmean = cell2mat(SIRData(3:end,10));
    data.GMMRsigma = cell2mat(SIRData(3:end,11));

    data.SRboundary = cell2mat(SIRData(3:end,13));
    data.SIboundary = cell2mat(SIRData(3:end,12));
    data.IRboundary = cell2mat(SIRData(3:end,14));
    data.Tecoff = cell2mat(SIRData(3:end,16));
    data.TecoffCIlow = cell2mat(SIRData(3:end,17));
    data.TecoffCIhi = cell2mat(SIRData(3:end,18));
    
    data.CLSIsticS = cell2mat(SIRData(3:end,19));
    data.CLSIsticR = cell2mat(SIRData(3:end,20));
    data.Xecoff = cell2mat(SIRData(3:end,21));
    data.XecoffCorrelation = cell2mat(SIRData(3:end,22));

    data.eucastBPS = cell2mat(SIRData(3:end,23));
    data.eucastBPR = cell2mat(SIRData(3:end,24));

    data.GMMweights = cell2mat(SIRData(3:end,25:27));

end