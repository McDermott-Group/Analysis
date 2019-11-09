function [nPhotons0, Qi0, Qc0] = plotMultifileExtractedQvsPower(n)
%Following function takes multiple files as inputs fits their resonances
%using Matt's function and outputs Qis and Qcs as function of photon number
nPhotons0 = zeros(1,n);  Qi0 = zeros(1,n);  Qc0 = zeros(1,n);

    for i=1:n
        data(i) = loadMeasurementData;
    end 
    for  i=1:n
        [nPhotons0(1,i), Qi0(1,i), Qc0(1,i), freq] = quarterWaveFitsLab1(data(i));
    end
        plti1 = figure;
        plti = axes('Parent',plti1);
        semilogx(nPhotons0,Qi0,'Marker','.','MarkerSize',20,'LineStyle','-.');
        set(plti,'XMinorTick','on','XScale','log');
        xlabel(['$\overline{n}','$'],'interpreter','latex','FontSize',16)
        ylabel('$Q_i$','interpreter','latex','FontSize',16)
        grid on
        
        pltc1 = figure;
        pltc = axes('Parent',pltc1);
        semilogx(nPhotons0,Qc0,'Marker','*','LineStyle','-.');
        set(pltc,'XMinorTick','on','XScale','log');
        xlabel(['$\overline{n}','$'],'interpreter','latex','FontSize',16)
        ylabel('$Q_c$','interpreter','latex','FontSize',16)
        grid on
end
    