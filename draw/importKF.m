function [Sow, X, Y, Z, B, L, H, E, N, U, Clock, Vx, Vy, Vz, VC, GS, BS] = importKF(filename, dataLines)
%IMPORTFILE 从文本文件中导入数据
%  [SOW, X, Y, Z, B, L, H, E, N, U, CLOCK, VX, VY, VZ, VC, GS, BS] =
%  IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  以列向量形式返回数据。
%
%  [SOW, X, Y, Z, B, L, H, E, N, U, CLOCK, VX, VY, VZ, VC, GS, BS] =
%  IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME 中的数据。对于不连续的行间隔，请将
%  DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  [Sow, X, Y, Z, B, L, H, E, N, U, Clock, Vx, Vy, Vz, VC, GS, BS] = importfile("D:\GitHub\Integrated-Navigation\draw\Static-KF.kf", [1, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2023-12-17 20:30:38 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [1, Inf];
end

%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 27);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = ["\t", ":"];

% 指定列名称和类型
opts.VariableNames = ["Var1", "Var2", "Sow", "Var4", "X", "Y", "Z", "Var8", "B", "L", "H", "Var12", "E", "N", "U", "Var16", "Clock", "Var18", "Vx", "Vy", "Vz", "VC", "Var23", "GS", "Var25", "BS", "Var27"];
opts.SelectedVariableNames = ["Sow", "X", "Y", "Z", "B", "L", "H", "E", "N", "U", "Clock", "Vx", "Vy", "Vz", "VC", "GS", "BS"];
opts.VariableTypes = ["string", "string", "double", "string", "double", "double", "double", "string", "double", "double", "double", "string", "double", "double", "double", "string", "double", "string", "double", "double", "double", "double", "string", "double", "string", "double", "string"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var8", "Var12", "Var16", "Var18", "Var23", "Var25", "Var27"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var8", "Var12", "Var16", "Var18", "Var23", "Var25", "Var27"], "EmptyFieldRule", "auto");

% 导入数据
tbl = readtable(filename, opts);

%% 转换为输出类型
Sow = tbl.Sow;
X = tbl.X;
Y = tbl.Y;
Z = tbl.Z;
B = tbl.B;
L = tbl.L;
H = tbl.H;
E = tbl.E;
N = tbl.N;
U = tbl.U;
Clock = tbl.Clock;
Vx = tbl.Vx;
Vy = tbl.Vy;
Vz = tbl.Vz;
VC = tbl.VC;
GS = tbl.GS;
BS = tbl.BS;
end