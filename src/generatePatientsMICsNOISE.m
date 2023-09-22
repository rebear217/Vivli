function [newdist,allnewPatientMICs] = generatePatientsMICsNOISE(dist,ldosages,p)

    M = length(ldosages);
    getdose = @(k)ldosages(k);

    patientN = sum(dist);
    allnewPatientMICs = zeros(1,patientN);

    pos = 1;
    for m = 1:M
        N = dist(m);
        if N > 0    
            R = 2*(rand(1,N)-0.5);
            if m == 1
                delta = 1*(R>1-p) + 1*(R>1-p*p);
            end
            if m == M
                delta = - 1*(R<-1+p) - 1*(R<-1+p*p);
            end
            if m > 1 && m < M
                delta = 1*(R>1-p) - 1*(R<-1+p);
            end
            if m > 2 && m < M-1
                delta = 1*(R>1-p*p) - 1*(R<-1+p*p);
            end

            newm = m + delta;
            newPatientMICs = getdose(newm);
            allnewPatientMICs(pos:pos+N-1) = newPatientMICs;
            pos = pos + N;
        end
    end

    newdist = histcounts(allnewPatientMICs,[ldosages 2*ldosages(end)]);

end