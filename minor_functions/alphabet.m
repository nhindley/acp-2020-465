
% alphabet.m, NPH 2017
% 
% USAGE:
% 
% alphabet(1) = 'a';
% alphabet(1,6) = {'a'}    {'b'}    {'c'}    {'d'}    {'e'}    {'f'}
% 

function OUT = alphabet(varargin)

alf = {
    'a'
    'b'
    'c'
    'd'
    'e'
    'f'
    'g'
    'h'
    'i'
    'j'
    'k'
    'l'
    'm'
    'n'
    'o'
    'p'
    'q'
    'r'
    's'
    't'
    'u'
    'v'
    'w'
    'x'
    'y'
    'z'};

% One input, output number/letter of alphabet that letter/number is
if nargin == 1
    IN = varargin{:};
    switch isnumeric(IN)
        case 1
            OUT = alf{IN};
        case 0
            OUT = find(strcmpi(alf,IN));
    end
end

% Two inputs, cell output of letters specified by the inputs
if nargin == 2
   st = varargin{1};
   ed = varargin{2};
   switch isnumeric(st)
       case 1
           OUT = {alf{st:ed}};
       case 0
           OUT = {alf{find(strcmpi(alf,st)):find(strcmpi(alf,ed))}};
   end
end

% should really have made this array safe, would be better input.
   
   
   
end
    