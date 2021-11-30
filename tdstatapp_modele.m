% Travaux dirig�s de Statistique Appliqu�e sous Matlab
% Script de d�part
clear all       % Efface les variables en m�moire
clc             % Efface le contenu de la fen�tre de commande
format compact  % Supprime les sauts de ligne dans la fen�tre de commande
cas = 'B';      % Exercice choisi (cf. le switch sur sa valeur)
ALPHA = 0.05;   % Risque d'erreur de premi�re esp�ce
chemin = '/Network/Servers/helvede.info.espci.fr/home/promo139/tpetrill/Documents/';    % Chemin pour les fichiers de donn�es, par exemple : '' en local, '/home/esa/irivals/MATLAB2A/' en salle info

close all
switch cas
    case 'A'
        display('Distributions de tailles (Galton)');
        noms = {'Pères', 'Méres', 'Couples'};
        nom_fichier = sprintf('%sdataA.txt', chemin);
        x = load(nom_fichier);
        x(:, 3) = (x(:, 1) + 1.08*x(:,2))/2
        for i = 1:3
            fprintf('\n*** %s ***\n', noms{i});
            %figure()
            tailles = x(:,i);
            
            % Tracer l'histogramme
            subplot(3,3,i)
            hist(tailles)
            title(sprintf('Histogram %s ', noms{i}))
            xlabel('Heights')
            ylabel('Effectives')
            
            % Tracer l'histogramme cumul�, etc.
            subplot(3,3,i+3)
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
            subplot(3,3,i+6)
            [h, p, stats] = chi2gof(tailles)
            qqplot(tailles);
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
        groupes = [zeros(10,1); ones(10,1)];
        subplot(2,2,1)
        boxplot(TRIOBPnorm, groupes, 'labels', {'TS21+', 'TS21-'})
        title('Boxplot TRIOBP')
        subplot(2,2,2)
        boxplot(MYLIPnorm, groupes, 'labels', {'TS21+', 'TS21-'})
        title('Boxplot MYLIP')        
        
        % Analyse diff�rentielle (comparaison des TS21+ et TS21-)
        noms = {'TRIOBP', 'MYLIP'};
        for i=1:2
            gene = data(:,i);
            
            %Ici on vais faire tous les tests
            geneplus = gene(1:10);
            genemoins = gene(11:20);
            
            [Hp, pValuep, SWstatisticp] = swtest(geneplus);
            [Hm, pValuem, SWstatisticm] = swtest(genemoins);

            if (Hp==0) && (Hm == 0)
                fprintf('\n*** G�ne %s ***\n', noms{i});
                fprintf('\n******* Test de la normalite *******\n')
               [hi,pi, statsi] = vartest2(geneplus, genemoins);
               if (hi == 0)
                   [hi,pi, statsi] = ttest2(geneplus, genemoins);
                   if (hi == 0)
                        fprintf('Egalite des esperances\n')
                   else
                        fprintf('Pas de egalite entre les esperances\n')
                   end
               else
                    [hi,pi, statsi] = ttest2(geneplus, genemoins, 'Vartype','unequal')
               end
            else
               fprintf('\n*** G�ne %s ***\n', noms{i});
               fprintf('\n******* Test nom parametrique *******\n')
               [pi,hi,statsi] = ranksum(geneplus, genemoins);
               if (hi == 0)
                   fprintf('Confirme le plus comme Gaussienne\n')
               else
                   fprintf('Non Gaussienne\n')
               end
            end
            
        end
end

