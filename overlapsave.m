function y = overlapsave(h, x, N)
% OVERLAPSAVE  实现 Overlap-Save (OLS) 快速卷积
%
% y = overlapsave(h, x, N)
%
% h : FIR 滤波器（长度 Lh）
% x : 输入信号（长度 Nx）
% N : FFT 点数，必须满足 N >= Lh
%
% 输出 y : 与 x 等长（Nx），为 OLS 线性卷积结果

h = h(:);
x = x(:);

Lh = length(h);          % 滤波器长度
Nx = length(x);

if N < Lh
    error('FFT 点数 N 必须 >= length(h)');
end

L = N - Lh + 1;          % 每块有效输出长度

% 输入前补 Lh-1 个零 (OLS 关键)
x_pad = [zeros(Lh-1,1); x];

% 滤波器频域
H = fft(h, N);

% 输出预分配
y = zeros(Nx, 1);

% 块处理
idx = 1;                 % 在 x_pad 中的读取位置
outpos = 1;              % 写输出的位置

while idx + N - 1 <= length(x_pad)
    
    % 取一块长度 N
    x_blk = x_pad(idx : idx+N-1);

    % FFT 卷积
    Y_blk = ifft( fft(x_blk) .* H );

    % 丢弃前 Lh-1 点（OLS 的 save）
    y_valid = Y_blk(Lh:end);

    % 写入输出
    y(outpos : min(outpos+L-1, Nx)) = ...
        y_valid(1 : min(L, Nx-outpos+1));

    % 下一块
    idx = idx + L;
    outpos = outpos + L;
end
end
