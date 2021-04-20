%load baseline_whisking4
group = ['abes'];   % abes or darb  
subgroup = ['abes'];  % abes or dar or db
animals = [1:20];    % abes 20 animals, dar 22, db 21
PODtrials = [0, 7, 14, 21, 28, 35, 42, 49, 56, 77, 92, 105];    % abes
%PODtrials = [0, 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,28,36,62];  % darb
% NB:  0 represents 'pre'
% pre=0 or POD7...14,21,28,35,42,49,56,77,92...POD105  (POD = post-operative
% day) for abes
% pre or POD 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,28,36,62 for darb

%PROBLEM:  abes.abes18POD14 .... chatter?
func1 = ['Blink'];  % Blink or Whisk
func2 = ['Whisk'];
side1 = ['L'];  % R or L
side2 = ['L'];
s1 = [func1 '-' side1];
s2 = [func2 '-' side2];

L = length(animals);
nT = length(PODtrials);
dynVAF = zeros(L,nT);
nFFT=200;
DF = 10;   % Set decimate factor

for n = 1:L
    animal = animals(n);
    trial = ['pre'];
    str1 = [group '.' subgroup, num2str(animal), '_', func1, trial, side1];
    str2 = [group '.' subgroup, num2str(animal), '_', func2, trial, side2];
    f1=decimate(eval(str1),DF);
    f2=decimate(eval(str2),DF);
    f1 = f1(1:6000);
    f2 = f2(1:6000);
    W = nldat(cat(2,f1,f2),'domainIncr',.01, 'chanNames',{s1 s2});
    W=ddt(W);
    W=W-mean(W);
    W=detrend(W); 
    W=ddt(W);
    W=smo(W,10);   % This really makes a BIG difference in VAF
    i=irf(W,'nLags',32,'nSides',2,'irfPseudoInvMode','auto')
    [R,V,yp]= nlid_resid(i,W,'plotFlag',false);
    dynVAF(n,1) = double(V);
    
    for m = 2:nT;
        trial = PODtrials(m);
        str1 = [group '.' subgroup, num2str(animal), '_', func1, 'POD', num2str(trial), side1];
        if isempty(eval(str1)) == 1
            dynVAF(n,m) == [];
        else
        str2 = [group '.' subgroup, num2str(animal), '_', func2, 'POD', num2str(trial), side2];
        f1=decimate(eval(str1),DF);
        f2=decimate(eval(str2),DF);
        f1 = f1(1:6000);
        f2 = f2(1:6000);
        W = nldat(cat(2,f1,f2),'domainIncr',.01, 'chanNames',{s1 s2});
        W=ddt(W);
        W=W-mean(W);
        W=detrend(W); 
        W=ddt(W);
        W=smo(W,10);   % This really makes a BIG difference in VAF
        i=irf(W,'nLags',32,'nSides',2,'irfPseudoInvMode','auto')
        [R,V,yp]= nlid_resid(i,W,'plotFlag',false);
        dynVAF(n,m) = double(V);
        end
    end
end


figure
stem(PODtrials, dynVAF', 'LineStyle', 'none');
hold on;
meanVAF = mean(dynVAF,1);
plot(PODtrials,meanVAF, '--r', 'LineWidth', 2);
hold off;
title(['Dynamic VAF: ', func1, '-', side1, ' vs. ', func2,'-',side2])
xlabel('Post operative day')
ylabel('VAF%')



