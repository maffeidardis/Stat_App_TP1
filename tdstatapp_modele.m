% Travaux dirig�s de Statistique Appliqu�e sous Matlab
% Script de d�part
clear all       % Efface les variables en m�moire
clc             % Efface le contenu de la fen�tre de commande
format compact  % Supprime les sauts de ligne dans la fen�tre de commande
cas = 'A';      % Exercice choisi (cf. le switch sur sa valeur)
ALPHA = 0.05;   % Risque d'erreur de premi�re esp�ce
chemin = '/Network/Servers/helvede.info.espci.fr/home/promo139/tpetrill/Documents/';    % Chemin pour les fichiers de donn�es, par exemple : '' en local, '/home/esa/irivals/MATLAB2A/' en salle info

close all
switch cas
    case 'A'
        display('Distributions de tailles (Galton)');
        noms = {'Pères', 'Méres', 'Couples'};
        nom_fichier = sprintf('%sdataA.txt', chemin);
        x = load(nom_fichier);
        for i = 1:2
            fprintf('\n*** %s ***\n', noms{i});
            %figure()
            tailles = x(:,i);
            
            % Tracer l'histogramme
            subplot(2,2,i)
            hist(tailles)
            title(sprintf('Histogram %s ', noms{i}))
            xlabel('Heights')
            ylabel('Effectives')
            
            % Tracer l'histogramme cumul�, etc.
            %subplot(2,2,i+2)
            [eff, xb] = hist(tailles);
            eff_sum = cumsum(eff);
            bar(xb, eff_sum);
            title(sprintf('Cumulative des effectives %s ', noms{i}))
            xlabel('Heights')
            ylabel('Cumulative Effectif')
            
            %Calcule des estimations biaises
            esperance = mean(tailles);
            variance = var(tailles);
            ecart_type = sqrt(variance);
            
            %Calcule des estimations non biaises
            variance_nb = var(tailles, 1);
            
            %Coefficients
            coeff_ass = skewness(tailles, 0);
            coeff_ass_nb = skewness(tailles, 1);
            coeff_aplatissement = kurtosis(tailles, 0);
            coeff_aplatissement_nb = kurtosis(tailles, 1);
            
            %Teste de caractère gaussien
            subplot(2,2,i+2)
            [h, p, stats] = chi2gof(tailles)
            
        end
        
    otherwise
        display('PCRq sur des donn�es de trisomie 21');
        nom_fichier = sprintf('%sdataB.txt', chemin);
        x = load(nom_fichier);

        % Moyenne des triplicats (i.e. moyennes des colonnes)
        HPRTmoy = mean(x(:,1:3),2);
        TRIOBPmoy = mean(x(:,4:6),2);
        MYLIPmoy  = mean(x(:,7:9),2);
        
        % Normaliser par rapport � HPRT (soustraction des ct)
        MYLIPnorm = MYLIPmoy-HPRTmoy;
        TRIOBPnorm = TRIOBPmoy-HPRTmoy;
        data = [TRIOBPnorm MYLIPnorm];
        
        % Tracer les bo�tes � moustaches
        
        
        % Analyse diff�rentielle (comparaison des TS21+ et TS21-)
        noms = {'TRIOBP', 'MYLIP'};
        for i=1:2
            fprintf('\n*** G�ne %s ***\n', noms{i});
            gene = data(:,i);
        end
end

