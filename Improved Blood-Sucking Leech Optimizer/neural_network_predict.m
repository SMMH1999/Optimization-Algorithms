%_________________________________________________________________________________
%  Ameliorated Golden jackal optimization (AGJO) with enhanced movement and multi-angle position updating strategy for solving engineering problems                                                                
%                                                                                                     
%  Developed in MATLAB R2022a                                                                  
%                                                                                                     
%  programming: Jianfu Bai                                                          
%                                                                                                     
%  e-Mail: Jianfu.Bai@UGent.be, magd.abdelwahab@ugent.be                                                               
%  Soete Laboratory, Department of Electrical Energy, Metals, Mechanical Constructions, and Systems, 
%  Faculty of Engineering and Architecture, Ghent University, Belgium                                                           
%                                                                                                                                                                                                                                                              
%  paper: Jianfu Bai, Samir Khatir, Laith Abualigah, Magd Abdel Wahab, Ameliorated Golden jackal optimization (AGJO) with enhanced movement and multi-angle position updating strategy for solving engineering problems (2023).  
%____________________________________________________________________________________

function predicted_output = neural_network_predict(X_train, weights_input_hidden, weights_hidden_output, biases_hidden, biases_output)

weights_input_hidden=reshape(weights_input_hidden, size(biases_hidden,2), size(X_train,2));
weights_hidden_output=reshape(weights_hidden_output, size(biases_output,2), size(biases_hidden,2));
z1 = weights_input_hidden * X_train' + biases_hidden';
a1 = sigmoid(z1);
z2 = weights_hidden_output * a1 + biases_output';
predicted_output = z2;

