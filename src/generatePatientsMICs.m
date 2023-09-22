function allnewPatientMICs = generatePatientsMICs(dist,ldosages)

    M = length(ldosages);

    patientN = sum(dist);
    allnewPatientMICs = zeros(1,patientN);

    pos = 1;
    for m = 1:M
        N = dist(m);
        d = ldosages(m);
        if N > 0    
            allnewPatientMICs(pos:pos+N-1) = d;
            pos = pos + N;
        end
    end

end