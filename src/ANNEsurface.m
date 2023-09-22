function ANNEsurface(bestNets,CEflag)

grey = [1 1 1]*0.3;
greyb = [1 1 2]*0.3;
greyr = [2 1 1]*0.3;

N = 500;
sigma = 1;
mu = [-3 2];

txt = 'EUCAST ANNE';
if CEflag
    txt = 'CLSI ANNE';
end

[h0,bins] = syntheticData(N,0,sigma);

h0 = nize(h0);
h1 = nize(syntheticData(N,mu(1),sigma));
h2 = nize(syntheticData(N,mu(2),sigma));

f = @(t,s)nize(s*(t*h0+(1-t)*h1)+(1-s)*h2);
%f = @(t,s)nize(s*h1+h0+t*h2);

colr = @(t,s,g)[t*(1-s/4),g,s*(1-t/4)];

filename = ['./mat/surfaceData_CEflag',num2str(CEflag),'.mat'];
if exist(filename,'file')
    load(filename,'BPLs','BPRs','synthdata','s','t')
    [M,~] = size(BPLs);
else
    M = 20;
    synthdata = cell(M);
    BPLs = zeros(M);
    BPRs = zeros(M);
    s = zeros(1,M);
    t = zeros(1,M);
    
    for i = 1:M
        t(i) = (i-1)/(M-1);
        for j = 1:M
            s(j) = (j-1)/(M-1);
            synthdata{i,j} = f(t(i),s(j));
            ANNbreakpointDecisions = ANNdecision(synthdata{i,j}',bestNets,[],CEflag);
            BPLs(i,j) = ANNbreakpointDecisions.bpLeftFinal;
            BPRs(i,j) = ANNbreakpointDecisions.bpRightFinal;
        end
        disp(num2str(t(i)));
    end
    save(filename,'BPLs','BPRs','synthdata','s','t')
end

[T,S] = meshgrid(t,s);
figure(1);
set(1,'pos',[123         629        2227         657])

subplot(3,3,[2 5 8])
m2 = meshc(T,S,BPRs,'DisplayName',[txt,' R breakpoint'],'LineWidth',1);
hold on
m2(1).EdgeColor = greyr;
m2(2).EdgeColor = greyr;
m2(2).HandleVisibility = 'off';
legend('location','northeast')
xlabel('t')
ylabel('s')
zlabel('R breakpoint ($\log_2\mu g/mL$)','Interpreter','latex')

subplot(3,3,[1 4 7])
m1 = meshc(T,S,BPLs,'DisplayName',[txt,' S breakpoint'],'LineWidth',1);
hold on
m1(1).EdgeColor = greyb;
m1(2).EdgeColor = greyb;
m1(2).HandleVisibility = 'off';
legend('location','northeast')
xlabel('t')
ylabel('s')
zlabel('S breakpoint ($\log_2\mu g/mL$)','Interpreter','latex')

for j = 1:1:M
    i = j;

    s1 = s(j);
    t1 = t(i);
    H1 = synthdata{i,j};
    bpl1 = BPLs(i,j);
    bpr1 = BPRs(i,j);

    i = 21-j;
    s2 = s(j);
    t2 = t(i);
    H2 = synthdata{i,j};
    bpl2 = BPLs(i,j);
    bpr2 = BPRs(i,j);
    
    i = 10;
    s3 = s(j);
    t3 = t(i);
    H3 = synthdata{i,j};
    bpl3 = BPLs(i,j);
    bpr3 = BPRs(i,j);
    
    colr1 = colr(t1,s1,0.25);
    colr2 = colr(t2,s2,0.75);
    colr3 = colr(2*t3,s3,1)*0.8;

    subplot(3,3,3)
    plot(bins,H1,'color',colr1,'linewidth',1,'HandleVisibility','off')
    hold on
    if (j==M)
        plot(bpl1,interp1(bins,H1,bpl1),'s','color',colr1/2,'markersize',2,'DisplayName','ANNE S breakpoints')
        plot(bpr1,interp1(bins,H1,bpr1),'.','color',colr1/2,'markersize',12,'DisplayName','ANNE R breakpoints')
    else
        plot(bpl1,interp1(bins,H1,bpl1),'s','color',colr1/2,'markersize',2,'HandleVisibility','off')
        plot(bpr1,interp1(bins,H1,bpr1),'.','color',colr1/2,'markersize',12,'HandleVisibility','off')
    end

    subplot(3,3,6)
    plot(bins,H2,'color',colr2,'linewidth',1,'HandleVisibility','off')
    hold on
    if (j==M)
        plot(bpl2,interp1(bins,H2,bpl2),'s','color',colr2/2,'markersize',2,'DisplayName','ANNE S breakpoints')
        plot(bpr2,interp1(bins,H2,bpr2),'.','color',colr2/2,'markersize',12,'DisplayName','ANNE R breakpoints')
    else
        plot(bpl2,interp1(bins,H2,bpl2),'s','color',colr2/2,'markersize',2,'HandleVisibility','off')
        plot(bpr2,interp1(bins,H2,bpr2),'.','color',colr2/2,'markersize',12,'HandleVisibility','off')
    end

    subplot(3,3,9)
    plot(bins,H3,'color',colr3,'linewidth',1,'HandleVisibility','off')
    hold on
    if (j==M)
        plot(bpl3,interp1(bins,H3,bpl3),'s','color',colr3/2,'markersize',2,'DisplayName','ANNE S breakpoints')
        plot(bpr3,interp1(bins,H3,bpr3),'.','color',colr3/2,'markersize',12,'DisplayName','ANNE R breakpoints')
    else
        plot(bpl3,interp1(bins,H3,bpl3),'s','color',colr3/2,'markersize',2,'HandleVisibility','off')
        plot(bpr3,interp1(bins,H3,bpr3),'.','color',colr3/2,'markersize',12,'HandleVisibility','off')
    end

    subplot(3,3,[2 5 8])
    if ~((j==1)||(j==M))
        plot3(s1,t1,bpr1,'.','HandleVisibility','off','markersize',30,'color',colr1)
        plot3(s2,t2,bpr2,'.','HandleVisibility','off','markersize',30,'color',colr2)
        plot3(s3,t3,bpr3,'.','HandleVisibility','off','markersize',30,'color',colr3)
    else
        plot3(s1,t1,bpr1,'.','DisplayName','path endpoint','markersize',30,'color',colr1)
        plot3(s2,t2,bpr2,'.','DisplayName','path endpoint','markersize',30,'color',colr2)
        plot3(s3,t3,bpr3,'.','DisplayName','path endpoint','markersize',30,'color',colr3)
    end

    subplot(3,3,[1 4 7])
    if ~((j==1)||(j==M))
        plot3(s1,t1,bpl1,'.','HandleVisibility','off','markersize',30,'color',colr1)
        plot3(s2,t2,bpl2,'.','HandleVisibility','off','markersize',30,'color',colr2)
        plot3(s3,t3,bpl3,'.','HandleVisibility','off','markersize',30,'color',colr3)
    else
        plot3(s1,t1,bpl1,'.','DisplayName','path endpoint','markersize',30,'color',colr1)
        plot3(s2,t2,bpl2,'.','DisplayName','path endpoint','markersize',30,'color',colr2)
        plot3(s3,t3,bpl3,'.','DisplayName','path endpoint','markersize',30,'color',colr3)
    end

end

subplot(3,3,[1 4 7])
view([-343,13])
zlim([-1/2 5.5])
text(1,0,0,'cyan')
text(0,1,3.5,'yellow')
text(0,0.5,3.5,'lime')
text(1,0.5,1.5,'green')
text(0,0,3.5,'black')
text(1,1,1.5,'purple')

subplot(3,3,[2 5 8])
view([-343,13])
zlim([-1/2 5.5])
text(1,0,2.25,'cyan')
text(0,1,5.25,'yellow')
text(0,0.5,5.25,'lime')
text(1,0.5,3.5,'green')
text(0,0,5.25,'black')
text(1,1,3.25,'purple')

subplot(3,3,3)
grid on
ylabel('freq.')
xlabel('MIC ($\log_2\mu g/mL$)','Interpreter','latex')
xlim([-6 6])
text(-5.5,0.35,'black-purple')
legend('boxoff')

subplot(3,3,6)
grid on
ylabel('freq.')
xlabel('MIC ($\log_2\mu g/mL$)','Interpreter','latex')
xlim([-6 6])
text(-5.5,0.35,'yellow-cyan')
legend('boxoff')

subplot(3,3,9)
grid on
ylabel('freq.')
xlabel('MIC ($\log_2\mu g/mL$)','Interpreter','latex')
xlim([-6 6])
text(-5.5,0.35,'lime-green')
legend('boxoff')

end

function H = nize(h)
    H = h/sum(h);
end

