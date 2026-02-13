function rv = calculate_rv(t, Max_iter)
    threshold = 0.1 * Max_iter;
    if t < threshold
        % 前期有一半概率使 rv 为 0.1
        if rand() < 0.5
            rv = 0.1;
        else
            rv = 0.9;  % 选择 0.9 以保持平均值接近 0.7
        end
    else
        % 后期使用非线性衰减，如指数衰减
        scale = 0.9 - 0.1;  % 衰减的范围，从 0.9 降到 0.1
        decay_rate = (Max_iter - threshold);  % 衰减的速度调节
        rv = 0.1 + scale * exp(-3 * (t - threshold) / decay_rate);
    end
end