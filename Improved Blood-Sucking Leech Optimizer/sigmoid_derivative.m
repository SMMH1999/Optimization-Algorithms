% 定义sigmoid函数的导数
function y = sigmoid_derivative(x)
    y = sigmoid(x) .* (1 - sigmoid(x));
end