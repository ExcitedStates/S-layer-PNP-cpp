function str = num2sci(num, vargin)

if nargin < 2
    pow10 = floor(log10(num));
    powsci = floor(pow10/3);
else
    powsci = vargin(1);
end

switch powsci
    case 4
        c = 'T';
        fac = 1e-12;
    case 3
        c = 'G';
        fac = 1e-9;
    case 2
        c = 'M';
        fac = 1e-6;
    case 1
        c = 'K';
        fac = 1e-3;
    case 0
        c = '';
        fac = 1;
    case -1
        c = 'm';
        fac = 1e3;
    case -2
        c = 'u';
        fac = 1e6;
    case -3
        c = 'n';
        fac = 1e9;
    case -4
        c = 'p';
        fac = 1e12;
    case -5
        c = 'f';
        fac = 1e15;
    case -6
        c = 'a';
        fac = 1e18;
    otherwise
        c = '';
        fac = 1;
end

str = [int0str( num*fac, 3) c];