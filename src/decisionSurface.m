function [Xs,Ys,Zs] = decisionSurface(entireSIRdata)

    data = extractSIRdataColumns(entireSIRdata);

    [Xs,Ys,Zs] = nanclean(data.GMMSmean,data.GMMRmean,data.CLSIsticS);
    [Xr,Yr,Zr] = nanclean(data.GMMSmean,data.GMMRmean,data.CLSIsticR);

    surffitS = fit([Xs,Ys],Zs,'poly44','normalize','on');
    surffitR = fit([Xr,Yr],Zr,'poly44','normalize','on');

    figure(1)
    plot3(Xs,Ys,Zs + 0.5*randn(size(Zs)),'.b','markersize',20)
    hold on
    plot3(Xr,Yr,Zr + 0.5*randn(size(Zr)),'.r','markersize',20)

    plot(surffitS,[Xs,Ys],Zs);
    colormap hsv
    plot(surffitR,[Xr,Yr],Zr);
    zlim([-9 9]);

    %{
    [Xs,Ys,Zs] = nanclean(data.SRboundary,data.IRboundary,data.CLSIsticS);
    [Xr,Yr,Zr] = nanclean(data.SRboundary,data.IRboundary,data.CLSIsticR);
    
    figure(2)
    plot3(Xs,Ys,Zs,'.b','markersize',20)
    hold on
    plot3(Xr,Yr,Zr,'.r','markersize',20)
    %}

    function [X,Y,Z] = nanclean(X,Y,Z)
        X(isnan(Y)) = 0;
        Y(isnan(Y)) = 0;

        Y(isnan(X)) = 0;
        X(isnan(X)) = 0;

        X = X(~isnan(Z));
        Y = Y(~isnan(Z));
        Z = Z(~isnan(Z));
 end
end