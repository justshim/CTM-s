function s = polyn2sym_mod(polyn)
% polyn2sympoly: convert a regression polynomial from polyfitn into its symbolic toolbox form
% usage: sp = polyn2sym(polyn)
%
% arguments: (input)
%  polyn - a structure as returned from polyfitn
%
% arguments: (output)
%  s - A symbolic toolbox object
%
% After conversion into a symbolic toolbox form, any symbolic operations are
% now possible on the polynomial.
%
% Requirement: The symbolic toolbox, as supplied by the MathWorks.
% http://www.mathworks.com/products/symbolic/functionlist.html
%
% See also: polyvaln, polyfit, polyval, sym
%
% Author: John D'Errico
% Release: 3.0
% Release date: 8/23/06

if exist('sym','file')~=2
    error 'Please obtain the symbolic toolbox from the MathWorks'
end

% initialize the returned argument as symbolic
%s = sym(0);

% Unpack the fields of polyn for use
Expon = polyn.ModelTerms;
Coef = polyn.Coefficients;
% how many terms?
nterms = length(Coef);

nvars = size(polyn.ModelTerms,2);

Varlist={};
for i = 1:nvars
    Varlist{i} = ['x(',num2str(i),')'];
end


% make the vars symbolic
% for i = 1:nvars
%     Varlist{i} = sym(Varlist{i});
% end

% build the polynomial
s = '@(x) ';
for i = 1:nterms
    %term = sym(Coef(i));
    term = num2str(Coef(i), '% +s');
    for j = 1:nvars
        if(Expon(i,j) ~= 0)
            term = strcat(term, '*', string(Varlist{j}), '^', num2str(Expon(i,j)));
        end
    end

    % accumulate into s
    s = strcat(s,term);
end

s = str2func(s);


