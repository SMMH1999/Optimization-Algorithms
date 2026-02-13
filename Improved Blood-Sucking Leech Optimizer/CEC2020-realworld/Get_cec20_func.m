function [lb,ub,dim,fobj] = Get_cec20_func(F)

switch F
    case 'F1'
        fobj = @Get_cec20_F1;
        dim = 9;                   % Dimension  
        lb = [0,0,0,0,1000,0,100,100,100];  % Lower boundary (vector)
        ub = [10,200,100,200,2000000,600,600,600,900];   % Uper boundary (vector)
        
    case 'F2'
        fobj = @Get_cec20_F2;
        dim = 11;                   % Dimension  
        lb = [10^4,10^4,10^4,0,0,0,100,100,100,100,100]; % Lower boundary (vector)
        ub = [0.819*10^6, 1.131*10^6, 2.05*10^6,0.05074,0.05074,0.05074,200,300,300,300,400];   % Uper boundary (vector)        
    case 'F3'
        fobj = @Get_cec20_F3;
        dim = 7;                   % Dimension  
        lb = [1000,0,2000,0,0,0,0];  % Lower boundary (vector)
        ub = [2000,100,4000,100,100,20,200];  % Uper boundary (vector)
        
    case 'F4'
        fobj = @Get_cec20_F4;
        dim = 6;                   % Dimension  
        lb = [0,0,0,0,1e-5,1e-5];  % Lower boundary (vector)
        ub = [1,1,1,1,16,16];
        
    case 'F5'
        fobj = @Get_cec20_F5;
        dim = 9;                   % Dimension  
        lb = -0*ones(1,dim);                
        ub= [100,200,100,100,100,100,200,100,200];
    case 'F6'
        fobj = @Get_cec20_F6;
        dim = 38;                   % Dimension  
        lb = 0*ones(1,dim);  % Lower boundary (vector)
        ub = [90,150,90,150,90,90,150,90,90,90,150,150,90,90,150,90,150,90,150,90,1,1.2,1,1,1,0.5,1,1,0.5,0.5,0.5,1.2,0.5,1.2,1.2,0.5,1.2,1.2];   % Uper boundary (vector)
        
    case 'F7'
        fobj = @Get_cec20_F7;
        dim = 48;                   % Dimension  
        lb = -0*ones(1,48); xmin7([24,26,28,31]) = 0.849999; % Lower boundary (vector)
        ub = 1*ones(1,dim); xmax7(4) = 140; xmax7([25,27,32,35,37,29]) = 30;xmax7([2,3,5,13,14,15]) = 90; xmax7([1,6,7,8,9,10,11,12,16,17,18,19,20]) = 35;   % Uper boundary (vector)        
    case 'F8'
        fobj = @Get_cec20_F8;
        dim = 2;                   % Dimension  
        lb = [0,-0.51];  % Lower boundary (vector)
        ub = [1.6,1.49];  % Uper boundary (vector)
        
    case 'F9'
        fobj = @Get_cec20_F9;
        dim = 3;                   % Dimension  
        lb = [0.5,0.5,-0.51];  % Lower boundary (vector)
        ub = [1.4,1.4,1.49];
        
    case 'F10'
        fobj = @Get_cec20_F10;
        dim = 3;                   % Dimension  
        lb = [0.2, -2.22554, -0.51];                
        ub= [1, -1, 1.49];      
    case 'F11'
        fobj = @Get_cec20_F11;
        dim = 7;                   % Dimension  
        lb = [0,0,0,0,-0.51,-0.51,0]; % Lower boundary (vector)
        ub = [20,20,10,10,1.49,1.49,40];  % Uper boundary (vector)
        
    case 'F12'
        fobj = @Get_cec20_F12;
        dim = 7;                   % Dimension  
        lb = [0,0,0,-0.51,-0.51,-0.51,-0.51];% Lower boundary (vector)
        ub = [100,100,100,1.49,1.49,1.49,1.49];   % Uper boundary (vector)        
    case 'F13'
        fobj = @Get_cec20_F13;
        dim = 5;                   % Dimension  
        lb = [27,27,27,77.51,32.51];  % Lower boundary (vector)
        ub = [45,45,45,102.49,45.49];  % Uper boundary (vector)
        
    case 'F14'
        fobj = @Get_cec20_F14;
        dim = 10;                   % Dimension  
        lb = [ 0.51,0.51,0.51,250,250,250,6,4,40,10];  % Lower boundary (vector)
        ub = [3.49,3.49,3.49,2500,2500,2500,20,16,700,450];
        
    case 'F15'
        fobj = @Get_cec20_F15;
        dim = 7;                   % Dimension  
        lb = [2.6, 0.7, 17, 7.3, 7.3, 2.9, 5];               
        ub= [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5];
    case 'F16'
        fobj = @Get_cec20_F16;
        dim = 14;                   % Dimension  
        lb = 0.001*ones(1,dim);  % Lower boundary (vector)
        ub = +5*ones(1,dim);   % Uper boundary (vector)
        
    case 'F17'
        fobj = @Get_cec20_F17;
        dim = 3;                   % Dimension  
        lb = [0.05,0.25,2.00]; % Lower boundary (vector)
        ub = [2,1.3,15.0];  % Uper boundary (vector)        
    case 'F18'
        fobj = @Get_cec20_F18;
        dim = 4;                   % Dimension  
        lb = [1,1,10,10];  % Lower boundary (vector)
        ub = [99,99,200,200];  % Uper boundary (vector)
        
    case 'F19'
        fobj = @Get_cec20_F19;
        dim = 4;                   % Dimension  
        lb = [0.125,0.1,0.1,0.1];  % Lower boundary (vector)
        ub = [2,10,10,2];
        
    case 'F20'
        fobj = @Get_cec20_F20;
        dim = 2;                   % Dimension  
        lb = 0*ones(1,dim);                
        ub= 1*ones(1,dim);    
    case 'F21'
        fobj = @Get_cec20_F21;
        dim = 5;                   % Dimension  
        lb = [60,90,1,0,2];  % Lower boundary (vector)
        ub = [80,110,3,1000,9];   % Uper boundary (vector)
        
    case 'F22'
        fobj = @Get_cec20_F22;
        dim = 9;                   % Dimension  
        lb = [16.51,13.51,13.51,16.51,13.51,47.51,0.51,0.51,0.51]; % Lower boundary (vector)
        ub = [96.49,54.49,51.49,46.49,51.49,124.49,3.49,6.49,6.49];   % Uper boundary (vector)        
    case 'F23'
        fobj = @Get_cec20_F23;
        dim = 5;                   % Dimension  
        lb = [0,0,0,0,0]; % Lower boundary (vector)
        ub =  [60,60,90,90,90];  % Uper boundary (vector)
        
    case 'F24'
        fobj = @Get_cec20_F24;
        dim = 7;                   % Dimension  
        lb = [10,10,100,0,10,100,1];  % Lower boundary (vector)
        ub = [150,150,200,50,150,300,3.14];
        
    case 'F25'
        fobj = @Get_cec20_F25;
        dim = 4;                   % Dimension  
        lb = [ 1, 1,  1e-6,1];                
        ub= [16, 16, 16*1e-6,16];
    case 'F26'
        fobj = @Get_cec20_F26;
        dim = 22;                   % Dimension  
        lb = [ 6.51.*ones(1,8), 0.51.*ones(1,14)];  % Lower boundary (vector)
        ub = [ 76.49.*ones(1,8), 4.49.*ones(1,4), 9.49.*ones(1,10)];   % Uper boundary (vector)
        
    case 'F27'
        fobj = @Get_cec20_F27;
        dim = 10;                   % Dimension  
        lb = 0.645e-4*ones(1,dim); % Lower boundary (vector)
        ub = 50e-4*ones(1,dim);   % Uper boundary (vector)        
    case 'F28'
        fobj = @Get_cec20_F28;
        dim = 10;                   % Dimension  
        lb = [125,10.5,4.51,0.515,0.515,0.4,0.6,0.3,0.02,0.6];  % Lower boundary (vector)
        ub =  [150,31.5,50.49,0.6,0.6,0.5,0.7,0.4,0.1,0.85];
        
    case 'F29'
        fobj = @Get_cec20_F29;
        dim = 4;                   % Dimension  
        lb = [20,1,20,0.1]; % Lower boundary (vector)
        ub =  [50,10,50,60];
        
    case 'F30'
        fobj = @Get_cec20_F30;
        dim = 3;                   % Dimension  
        lb = [0.51,0.6,0.51];               
        ub= [70.49,3,42.49];      
    case 'F31'
        fobj = @Get_cec20_F31;
        dim = 4;                   % Dimension  
        lb = 12.*ones(1,4); % Lower boundary (vector)
        ub = 60.*ones(1,4);  % Uper boundary (vector)
        
    case 'F32'
        fobj = @Get_cec20_F32;
        dim = 5;                   % Dimension  
        lb = [78,33,27,27,27];% Lower boundary (vector)
        ub = [102,45,45,45,45];   % Uper boundary (vector)        
    case 'F33'
        fobj = @Get_cec20_F33;
        dim = 30;                   % Dimension  
        lb = 0.001.*ones(1,dim);  % Lower boundary (vector)
        ub = ones(1,dim);  % Uper boundary (vector)
        
    case 'F34'
        fobj = @Get_cec20_F34;
        dim = 118;                   % Dimension  
        lb = -1*ones(1,dim); % Lower boundary (vector)
        ub = +1*ones(1,dim);
        
    case 'F35'
        fobj = @Get_cec20_F35;
        dim = 153;                   % Dimension  
        lb = -1*ones(1,dim);               
        ub= +1*ones(1,dim);
    case 'F36'
        fobj = @Get_cec20_F36;
        dim = 158;                   % Dimension  
        lb = -1*ones(1,dim);  % Lower boundary (vector)
        ub = +1*ones(1,dim);   % Uper boundary (vector)
        
    case 'F37'
        fobj = @Get_cec20_F37;
        dim = 126;                   % Dimension  
        lb = -1*ones(1,dim);xmin37(117:126) = 0;
        ub = +1*ones(1,dim);% Uper boundary (vector)        
    case 'F38'
        fobj = @Get_cec20_F38;
        dim = 126;                   % Dimension  
        lb = -1*ones(1,dim);xmin38(117:126) = 0;  % Lower boundary (vector)
        ub = +1*ones(1,dim);  % Uper boundary (vector)
        
    case 'F39'
        fobj = @Get_cec20_F39;
        dim = 126;                   % Dimension  
        lb = -1*ones(1,dim);xmin39(117:126) = 0;  % Lower boundary (vector)
        ub = +1*ones(1,dim);
        
    case 'F40'
        fobj = @Get_cec20_F40;
        dim = 76;                   % Dimension  
        lb = -1*ones(1,dim);xmin40(75:76) = 0;                
        ub= +1*ones(1,dim);xmax40(75:76) = 2;      
    case 'F41'
        fobj = @Get_cec20_F41;
        dim = 74;                   % Dimension  
        lb = -1*ones(1,dim);  % Lower boundary (vector)
        ub = +1*ones(1,dim);   % Uper boundary (vector)
        
    case 'F42'
        fobj = @Get_cec20_F42;
        dim = 86;                   % Dimension  
        lb = -1*ones(1,dim);xmin42(75:76) = 0;xmin42(77:86) = 0; % Lower boundary (vector)
        ub = +1*ones(1,dim);xmax42(75:76) = 2;xmax42(77:86) = 500;   % Uper boundary (vector)        
    case 'F43'
        fobj = @Get_cec20_F43;
        dim = 86;                   % Dimension  
        lb = -1*ones(1,dim);xmin43(75:76) = 0;xmin43(77:86) = 0;  % Lower boundary (vector)
        ub = +1*ones(1,dim);xmax43(75:76) = 2;xmax43(77:86) = 500;  % Uper boundary (vector)
        
    case 'F44'
        fobj = @Get_cec20_F44;
        dim = 30;                   % Dimension  
        lb = 40*ones(1,dim);  % Lower boundary (vector)
        ub = 1960*ones(1,dim);
        
    case 'F45'
        fobj = @Get_cec20_F45;
        dim = 25;                   % Dimension  
        lb = -0*ones(1,dim);               
        ub= +90*ones(1,dim);
    case 'F46'
        fobj = @Get_cec20_F46;
        dim = 25;                   % Dimension  
        lb = -0*ones(1,dim);  % Lower boundary (vector)
        ub = +90*ones(1,dim);
        
    case 'F47'
        fobj = @Get_cec20_F47;
        dim = 25;                   % Dimension  
        lb = -0*ones(1,dim); % Lower boundary (vector)
        ub = +90*ones(1,dim);
    case 'F48'
        fobj = @Get_cec20_F48;
        dim = 30;                   % Dimension  
        lb = -0*ones(1,dim);  % Lower boundary (vector)
        ub = +90*ones(1,dim);  % Uper boundary (vector)
        
    case 'F49'
        fobj = @Get_cec20_F49;
        dim = 30;                   % Dimension  
        lb = -0*ones(1,dim); % Lower boundary (vector)
        ub = +90*ones(1,dim);
        
    case 'F50'
        fobj = @Get_cec20_F50;
        dim = 30;                   % Dimension  
        lb = -0*ones(1,dim);           
        ub= +90*ones(1,dim); 
    case 'F51'
        fobj = @Get_cec20_F51;
        dim = 59;                   % Dimension  
        lb = 0.*ones(1,dim); % Lower boundary (vector)
        ub = 10.*ones(1,dim); % Uper boundary (vector)
        
    case 'F52'
        fobj = @Get_cec20_F52;
        dim = 59;                   % Dimension  
        lb = 0.*ones(1,dim);% Lower boundary (vector)
        ub = 10.*ones(1,dim);  % Uper boundary (vector)        
    case 'F53'
        fobj = @Get_cec20_F53;
        dim = 59;                   % Dimension  
        lb = 0.*ones(1,dim);  % Lower boundary (vector)
        ub = 10.*ones(1,dim);  % Uper boundary (vector)
        
    case 'F54'
        fobj = @Get_cec20_F54;
        dim = 59;                   % Dimension  
        lb = 0.*ones(1,dim);  % Lower boundary (vector)
        ub = 10.*ones(1,dim);
        
    case 'F55'
        fobj = @Get_cec20_F55;
        dim = 64;                   % Dimension  
        lb = 0.*ones(1,dim);               
        ub= 10.*ones(1,dim);
    case 'F56'
        fobj = @Get_cec20_F56;
        dim = 64;                   % Dimension  
        lb = 0.*ones(1,dim);  % Lower boundary (vector)
        ub = 10.*ones(1,dim);   % Uper boundary (vector)
        
    case 'F57'
        fobj = @Get_cec20_F57;
        dim = 64;                   % Dimension  
        lb = 0.*ones(1,dim); % Lower boundary (vector)
        ub = 10.*ones(1,dim);  % Uper boundary (vector)                
end

end

%% Industrial Chemical Processes		
function o = Get_cec20_F1(x)
    %% Heat Exchanger Network Design (case 1)
    [ps,D]=size(x);
    f = 35.*x(:,1).^0.6 + 35.*x(:,2).^0.6;

    h1 = 200.*x(:,1).*x(:,4)-x(:,3);
    h2 = 200.*x(:,2).*x(:,6)-x(:,5);
    h3 = x(:,3) - 10000.*(x(:,7)-100);
    h4 = x(:,5) - 10000.*(300-x(:,7));
    h5 = x(:,3) - 10000.*(600-x(:,8));
    h6 = x(:,5) - 10000.*(900-x(:,9));
    h7 = x(:,4).*log(abs(x(:,8)-100)+1e-8)-x(:,4).*log((600-x(:,7))+1e-8)-x(:,8)+x(:,7)+500;
    h8 = x(:,6).*log(abs(x(:,9)-x(:,7))+1e-8)-x(:,6).*log(600)-x(:,9)+x(:,7)+600;
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    H5 =abs(h5)-0.0001;
    H6 =abs(h6)-0.0001;
    H7 =abs(h7)-0.0001;
    H8 =abs(h8)-0.0001;    
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;
    panalty_5 = 10e100*(max(0,H5))^2;
    panalty_6 = 10e100*(max(0,H6))^2;
    panalty_7 = 10e100*(max(0,H7))^2;
    panalty_8 = 10e100*(max(0,H8))^2;

    o  = f + panalty_1+panalty_2+panalty_3+panalty_4+panalty_5+panalty_6+panalty_7+panalty_8;
end


function o = Get_cec20_F2(x)
    %% Heat Exchanger Network Design (case 2)
    [ps,D]=size(x);
    f = (x(:,1)./(120*x(:,4))).^0.6+(x(:,2)./(80*x(:,5))).^0.6+(x(:,3)./(40*x(:,6))).^0.6;

    h1 = x(:,1)-1e4.*(x(:,7)-100);
    h2 = x(:,2)-1e4.*(x(:,8)-x(7));
    h3 = x(:,3)-1e4.*(500-x(:,8));
    h4 = x(:,1)-1e4.*(300-x(:,9));
    h5 = x(:,2)-1e4.*(400-x(:,10));
    h6 = x(:,3)-1e4.*(600-x(:,11));
    h7 = x(:,4).*log(abs(x(:,9)-100)+1e-8)-x(:,4).*log(300-x(:,7)+1e-8)-x(:,9)-x(:,7)+400;
    h8 = x(:,5).*log(abs(x(:,10)-x(:,7))+1e-8)-x(:,5).*log(abs(400-x(:,8))+1e-8)-x(:,10)+x(:,7)-x(:,8)+400;
    h9 = x(:,6).*log(abs(x(:,11)-x(:,8))+1e-8)-x(:,6).*log(100)-x(:,11)+x(:,8)+100;
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    H5 =abs(h5)-0.0001;
    H6 =abs(h6)-0.0001;
    H7 =abs(h7)-0.0001;
    H8 =abs(h8)-0.0001; 
    H9 =abs(h9)-0.0001;
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;
    panalty_5 = 10e100*(max(0,H5))^2;
    panalty_6 = 10e100*(max(0,H6))^2;
    panalty_7 = 10e100*(max(0,H7))^2;
    panalty_8 = 10e100*(max(0,H8))^2;    
    panalty_9 = 10e100*(max(0,H9))^2;
    o  = f + panalty_1+panalty_2+panalty_3+panalty_4+panalty_5+panalty_6+panalty_7+panalty_8+panalty_9;    
end

function o = Get_cec20_F3(x)
    %% Optimal Operation of Alkylation Unit
    [ps,D]=size(x);
      f = -1.715.*x(:,1)-0.035.*x(:,1).*x(:,6)-4.0565.*x(:,3)-10.0.*x(:,2)+0.063.*x(:,3).*x(:,5);

      g1 = 0.0059553571.*x(:,6).^2.*x(:,1)+0.88392857.*x(:,3)-0.1175625.*x(:,6).*x(:,1)-x(:,1);
      g2 = 1.1088.*x(:,1)+0.1303533.*x(:,1).*x(:,6)-0.0066033.*x(:,1).*x(:,6).^2-x(:,3);
      g3 = 6.66173269.*x(:,6).^2+172.39878.*x(:,5)-56.596669.*x(:,4)-191.20592.*x(:,6)-10000;
      g4 = 1.08702.*x(:,6)+0.32175.*x(:,4)-0.03762.*x(:,6).^2-x(:,5)+56.85075;
      g5 = 0.006198.*x(:,7).*x(:,4).*x(:,3)+2462.3121.*x(:,2)-25.125634.*x(:,2).*x(:,4)-x(:,3).*x(:,4);
      g6 = 161.18996.*x(:,3).*x(:,4)+5000.0.*x(:,2).*x(:,4)-489510.0.*x(:,2)-x(:,3).*x(:,4).*x(:,7);
      g7 = 0.33.*x(:,7)-x(:,5)+44.333333;
      g8 = 0.022556.*x(:,5)-0.007595.*x(:,7)-1.0;
      g9 = 0.00061.*x(:,3)-0.0005.*x(:,1)-1.0;
      g10= 0.819672.*x(:,1)-x(:,3)+0.819672;
      g11= 24500.0.*x(:,2)-250.0.*x(:,2).*x(:,4)-x(:,3).*x(:,4);
      g12= 1020.4082.*x(:,4).*x(:,2)+1.2244898.*x(:,3).*x(:,4)-100000.*x(:,2);
      g13= 6.25.*x(:,1).*x(:,6)+6.25.*x(:,1)-7.625.*x(:,3)-100000;
      g14= 1.22.*x(:,3)-x(:,6).*x(:,1)-x(:,1)+1.0;
      panalty_1 = 10e100*(max(0,g1))^2;
      panalty_2 = 10e100*(max(0,g2))^2;
      panalty_3 = 10e100*(max(0,g3))^2;
      panalty_4 = 10e100*(max(0,g4))^2;
      panalty_5 = 10e100*(max(0,g5))^2;
      panalty_6 = 10e100*(max(0,g6))^2;
      panalty_7 = 10e100*(max(0,g7))^2;
      panalty_8 = 10e100*(max(0,g8))^2;
      panalty_9 = 10e100*(max(0,g9))^2;
      panalty_10 = 10e100*(max(0,g10))^2;
      panalty_11 = 10e100*(max(0,g11))^2;
      panalty_12 = 10e100*(max(0,g12))^2;
      panalty_13 = 10e100*(max(0,g13))^2;
      panalty_14 = 10e100*(max(0,g14))^2;
      o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14;       
end

function o = Get_cec20_F4(x)
    %% Reactor Network Design (RND)
    [ps,D]=size(x);
    k1 = 0.09755988;
    k2 = 0.99.*k1;
    k3 = 0.0391908;
    k4 = 0.9.*k3;
    f = -x(:,4);
    h1 = x(:,1)+k1.*x(:,2).*x(:,5)-1;
    h2 = x(:,2)-x(:,1)+k2.*x(:,2).*x(:,6);
    h3 = x(:,3)+x(:,1)+k3.*x(:,3).*x(:,5)-1;
    h4 = x(:,4)-x(:,3)+x(:,2)-x(:,1)+k4.*x(:,4).*x(:,6);
    g1 = x(:,5).^0.5+x(:,6).^0.5-4;
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;    
    panalty_5 = 10e100*(max(0,g1))^2;
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5;
end

function o = Get_cec20_F5(x)
    %% Haverly's Pooling Problem
    [ps,D]=size(x);
    f = -(9.*x(:,1)+15.*x(:,2)-6.*x(:,3)-16.*x(:,4)-10.*(x(:,5)+x(:,6)));
    g1 = x(:,9).*x(:,7)+2.*x(:,5)-2.5.*x(:,1);
    g2 = x(:,9).*x(:,8)+2.*x(:,6)-1.5.*x(:,2);
    h1 = x(:,7)+x(:,8)-x(:,3)-x(:,4);
    h2 = x(:,1)-x(:,7)-x(:,5);
    h3 = x(:,2)-x(:,8)-x(:,6);
    h4 = x(:,9).*x(:,7)+x(:,9).*x(:,8)-3.*x(:,3)-x(:,4);
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    panalty_3 = 10e100*(max(0,H1))^2;
    panalty_4 = 10e100*(max(0,H2))^2;
    panalty_5 = 10e100*(max(0,H3))^2;
    panalty_6 = 10e100*(max(0,H4))^2;
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6;
end

function o = Get_cec20_F6(x)
    %% Blending-Pooling-Separation problem
    [ps,D]=size(x);
    f = 0.9979+0.00432.*x(:,5)+0.01517.*x(:,13);

    h1 = x(:,1)+x(:,2)+x(:,3)+x(:,4)-300;
    h2 = x(:,6)-x(:,7)-x(:,8);
    h3 = x(:,9)-x(:,10)-x(:,11)-x(:,12);
    h4 = x(:,14)-x(:,15)-x(:,16)-x(:,17);
    h5 = x(:,18)-x(:,19)-x(:,20);
    h6 = x(:,5).*x(:,21)-x(:,6).*x(:,22)-x(:,9).*x(:,23);
    h7 = x(:,5).*x(:,24)-x(:,6).*x(:,25)-x(:,9).*x(:,26);
    h8 = x(:,5).*x(:,27)-x(:,6).*x(:,28)-x(:,9).*x(:,29);
    h9 = x(:,13).*x(:,30)-x(:,14).*x(:,31)-x(:,18).*x(:,32);
    h10 = x(:,13).*x(:,33)-x(:,14).*x(:,34)-x(:,18).*x(:,35);
    h11 = x(:,13).*x(:,36)-x(:,14).*x(:,37)-x(:,18).*x(:,38);
    h12 = 1/3.*x(:,1)+x(:,15).*x(:,31)-x(:,5).*x(:,21);
    h13 = 1/3.*x(:,1)+x(:,15).*x(:,34)-x(:,5).*x(:,24);
    h14 = 1/3.*x(:,1)+x(:,15).*x(:,37)-x(:,5).*x(:,27);
    h15 = 1/3.*x(:,2)+x(:,10).*x(:,23)-x(:,13).*x(:,30);
    h16 = 1/3.*x(:,2)+x(:,10).*x(:,26)-x(:,13).*x(:,33);
    h17 = 1/3.*x(:,2)+x(:,10).*x(:,29)-x(:,13).*x(:,36);
    h18 = 1/3.*x(:,3)+x(:,7).*x(:,22)+x(:,11).*x(:,23)+x(:,16).*x(:,31)+x(:,19).*x(:,32)-30;
    h19 = 1/3.*x(:,3)+x(:,7).*x(:,25)+x(:,11).*x(:,26)+x(:,16).*x(:,34)+x(:,19).*x(:,35)-50;
    h20 = 1/3.*x(:,3)+x(:,7).*x(:,28)+x(:,11).*x(:,29)+x(:,16).*x(:,37)+x(:,19).*x(:,38)-30;
    h21 = x(:,21)+x(:,24)+x(:,27)-1;
    h22 = x(:,22)+x(:,25)+x(:,28)-1;
    h23 = x(:,23)+x(:,26)+x(:,29)-1;
    h24 = x(:,30)+x(:,33)+x(:,36)-1;
    h25 = x(:,31)+x(:,34)+x(:,37)-1;
    h26 = x(:,32)+x(:,35)+x(:,38)-1;
    h27 = x(:,25);
    h28 = x(:,28);
    h29 = x(:,23);
    h30 = x(:,37);
    h31 = x(:,32);
    h32 = x(:,35);
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    H5 =abs(h5)-0.0001;
    H6 =abs(h6)-0.0001;
    H7 =abs(h7)-0.0001;
    H8 =abs(h8)-0.0001; 
    H9 =abs(h9)-0.0001;
    H10 =abs(h10)-0.0001;
    H11 =abs(h11)-0.0001;
    H12 =abs(h12)-0.0001;
    H13 =abs(h13)-0.0001;
    H14 =abs(h14)-0.0001;
    H15 =abs(h15)-0.0001;
    H16 =abs(h16)-0.0001;
    H17 =abs(h17)-0.0001;
    H18 =abs(h18)-0.0001; 
    H19 =abs(h19)-0.0001;
    H20 =abs(h20)-0.0001; 
    H21 =abs(h21)-0.0001;
    H22 =abs(h22)-0.0001;
    H23 =abs(h23)-0.0001;
    H24 =abs(h24)-0.0001;
    H25 =abs(h25)-0.0001;
    H26 =abs(h26)-0.0001;
    H27 =abs(h27)-0.0001;
    H28 =abs(h28)-0.0001; 
    H29 =abs(h29)-0.0001;
    H30 =abs(h30)-0.0001;
    H31 =abs(h31)-0.0001;
    H32 =abs(h32)-0.0001;
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;
    panalty_5 = 10e100*(max(0,H5))^2;
    panalty_6 = 10e100*(max(0,H6))^2;
    panalty_7 = 10e100*(max(0,H7))^2;
    panalty_8 = 10e100*(max(0,H8))^2;    
    panalty_9 = 10e100*(max(0,H9))^2;  
    panalty_10 = 10e100*(max(0,H10))^2; 
    panalty_11 = 10e100*(max(0,H11))^2;
    panalty_12 = 10e100*(max(0,H12))^2;
    panalty_13 = 10e100*(max(0,H13))^2;
    panalty_14 = 10e100*(max(0,H14))^2;
    panalty_15 = 10e100*(max(0,H15))^2;
    panalty_16 = 10e100*(max(0,H16))^2;
    panalty_17 = 10e100*(max(0,H17))^2;
    panalty_18 = 10e100*(max(0,H18))^2;    
    panalty_19 = 10e100*(max(0,H19))^2;  
    panalty_20 = 10e100*(max(0,H20))^2; 
    panalty_21 = 10e100*(max(0,H21))^2;
    panalty_22 = 10e100*(max(0,H22))^2;
    panalty_23 = 10e100*(max(0,H23))^2;
    panalty_24 = 10e100*(max(0,H24))^2;
    panalty_25 = 10e100*(max(0,H25))^2;
    panalty_26 = 10e100*(max(0,H26))^2;
    panalty_27 = 10e100*(max(0,H27))^2;
    panalty_28 = 10e100*(max(0,H28))^2;    
    panalty_29 = 10e100*(max(0,H29))^2;  
    panalty_30 = 10e100*(max(0,H30))^2; 
    panalty_31 = 10e100*(max(0,H31))^2;  
    panalty_32 = 10e100*(max(0,H32))^2;   
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15 + panalty_16 + panalty_17 + panalty_18 + panalty_19...
         + panalty_20 + panalty_21 + panalty_22 + panalty_23 + panalty_24 + panalty_25 + panalty_26 + panalty_27 + panalty_28 + panalty_29...
         + panalty_30  + panalty_31 + panalty_32;      
end

function o = Get_cec20_F7(x)
    %% Propane, Isobutane, n-Butane Nonsharp Separation
    [ps,D]=size(x);
    c = [ 0.23947, 0.75835; -0.0139904, -0.0661588; 0.0093514, 0.0338147;.....
        0.0077308, 0.0373349; -0.0005719, 0.0016371;0.0042656, 0.0288996];
    f = c(1,1)+(c(2,1)+c(3,1).*x(:,24)+c(4,1).*x(:,28)+c(5,1).*x(:,33)+c(6,1).*x(:,34)).*x(:,5)....
        +c(1,2)+(c(2,2)+c(3,2).*x(:,26)+c(4,2).*x(:,31)+c(5,2).*x(:,38)+c(6,2).*x(:,39)).*x(:,13);

    h1 = x(:,1)+x(:,2)+x(:,3)+x(:,4)-300;
    h2 = x(:,6)-x(:,7)-x(:,8);
    h3 = x(:,9)-x(:,10)-x(:,11)-x(:,12);
    h4 = x(:,14)-x(:,15)-x(:,16)-x(:,17);
    h5 = x(:,18)-x(:,19)-x(:,20);
    h6 = x(:,6).*x(:,21)-x(:,24).*x(:,25);
    h7 = x(:,14).*x(:,22)-x(:,26).*x(:,27);
    h8 = x(:,9).*x(:,23)-x(:,28).*x(:,29);
    h9 = x(:,18).*x(:,30)-x(:,31).*x(:,32);
    h10 = x(:,25)-x(:,5).*x(:,33);
    h11 = x(:,29)-x(:,5).*x(:,34);
    h12 = x(:,35)-x(:,5).*x(:,36);
    h13 = x(:,37)-x(:,13).*x(:,38);
    h14 = x(:,27)-x(:,13).*x(:,39);
    h15 = x(:,32)-x(:,13).*x(:,40);
    h16 = x(:,25)-x(:,6).*x(:,21)-x(:,9).*x(:,41);
    h17 = x(:,29)-x(:,6).*x(:,42)-x(:,9).*x(:,23);
    h18 = x(:,35)-x(:,6).*x(:,43)-x(:,9).*x(:,44);
    h19 = x(:,37)-x(:,14).*x(:,45)-x(:,18).*x(:,46);
    h20 = x(:,27)-x(:,14).*x(:,22)-x(:,18).*x(:,47);
    h21 = x(:,32)-x(:,14).*x(:,48)-x(:,18).*x(:,30);
    h22 = 1/3*x(:,1)+x(:,15).*x(:,45)-x(:,25);
    h23 = 1/3*x(:,1)+x(:,15).*x(:,22)-x(:,29);
    h24 = 1/3*x(:,1)+x(:,15).*x(:,48)-x(:,35);
    h25 = 1/3*x(:,2)+x(:,10).*x(:,41)-x(:,37);
    h26 = 1/3*x(:,2)+x(:,10).*x(:,23)-x(:,27);
    h27 = 1/3*x(:,2)+x(:,10).*x(:,44)-x(:,32);
    h28 = x(:,33)+x(:,34)+x(:,36)-1;
    h29 = x(:,21)+x(:,42)+x(:,43)-1;
    h30 = x(:,41)+x(:,23)+x(:,44)-1;
    h31 = x(:,38)+x(:,39)+x(:,40)-1;
    h32 = x(:,45)+x(:,22)+x(:,48)-1;
    h33 = x(:,46)+x(:,47)+x(:,30)-1;
    h34 = x(:,43);
    h35 = x(:,46);
    h36 = 1/3*x(:,3)+x(:,7).*x(:,21)+x(:,11).*x(:,41)+x(:,16).*x(:,45)+x(:,19).*x(:,46)-30;
    h37 = 1/3*x(:,3)+x(:,7).*x(:,42)+x(:,11).*x(:,23)+x(:,16).*x(:,22)+x(:,19).*x(:,47)-50;
    h38 = 1/3*x(:,3)+x(:,7).*x(:,43)+x(:,11).*x(:,44)+x(:,16).*x(:,48)+x(:,19).*x(:,30)-30;
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    H5 =abs(h5)-0.0001;
    H6 =abs(h6)-0.0001;
    H7 =abs(h7)-0.0001;
    H8 =abs(h8)-0.0001; 
    H9 =abs(h9)-0.0001;
    H10 =abs(h10)-0.0001;
    H11 =abs(h11)-0.0001;
    H12 =abs(h12)-0.0001;
    H13 =abs(h13)-0.0001;
    H14 =abs(h14)-0.0001;
    H15 =abs(h15)-0.0001;
    H16 =abs(h16)-0.0001;
    H17 =abs(h17)-0.0001;
    H18 =abs(h18)-0.0001; 
    H19 =abs(h19)-0.0001;
    H20 =abs(h20)-0.0001; 
    H21 =abs(h21)-0.0001;
    H22 =abs(h22)-0.0001;
    H23 =abs(h23)-0.0001;
    H24 =abs(h24)-0.0001;
    H25 =abs(h25)-0.0001;
    H26 =abs(h26)-0.0001;
    H27 =abs(h27)-0.0001;
    H28 =abs(h28)-0.0001; 
    H29 =abs(h29)-0.0001;
    H30 =abs(h30)-0.0001;
    H31 =abs(h31)-0.0001;
    H32 =abs(h32)-0.0001;
    H33 =abs(h33)-0.0001;
    H34 =abs(h34)-0.0001;
    H35 =abs(h35)-0.0001;
    H36 =abs(h36)-0.0001;
    H37 =abs(h37)-0.0001;
    H38 =abs(h38)-0.0001;    
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;
    panalty_5 = 10e100*(max(0,H5))^2;
    panalty_6 = 10e100*(max(0,H6))^2;
    panalty_7 = 10e100*(max(0,H7))^2;
    panalty_8 = 10e100*(max(0,H8))^2;    
    panalty_9 = 10e100*(max(0,H9))^2;  
    panalty_10 = 10e100*(max(0,H10))^2; 
    panalty_11 = 10e100*(max(0,H11))^2;
    panalty_12 = 10e100*(max(0,H12))^2;
    panalty_13 = 10e100*(max(0,H13))^2;
    panalty_14 = 10e100*(max(0,H14))^2;
    panalty_15 = 10e100*(max(0,H15))^2;
    panalty_16 = 10e100*(max(0,H16))^2;
    panalty_17 = 10e100*(max(0,H17))^2;
    panalty_18 = 10e100*(max(0,H18))^2;    
    panalty_19 = 10e100*(max(0,H19))^2;  
    panalty_20 = 10e100*(max(0,H20))^2; 
    panalty_21 = 10e100*(max(0,H21))^2;
    panalty_22 = 10e100*(max(0,H22))^2;
    panalty_23 = 10e100*(max(0,H23))^2;
    panalty_24 = 10e100*(max(0,H24))^2;
    panalty_25 = 10e100*(max(0,H25))^2;
    panalty_26 = 10e100*(max(0,H26))^2;
    panalty_27 = 10e100*(max(0,H27))^2;
    panalty_28 = 10e100*(max(0,H28))^2;    
    panalty_29 = 10e100*(max(0,H29))^2;  
    panalty_30 = 10e100*(max(0,H30))^2; 
    panalty_31 = 10e100*(max(0,H31))^2;  
    panalty_32 = 10e100*(max(0,H32))^2; 
    panalty_33 = 10e100*(max(0,H33))^2; 
    panalty_34 = 10e100*(max(0,H34))^2;  
    panalty_35 = 10e100*(max(0,H35))^2; 
    panalty_36 = 10e100*(max(0,H36))^2; 
    panalty_37 = 10e100*(max(0,H37))^2;  
    panalty_38 = 10e100*(max(0,H38))^2;     
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15 + panalty_16 + panalty_17 + panalty_18 + panalty_19...
         + panalty_20 + panalty_21 + panalty_22 + panalty_23 + panalty_24 + panalty_25 + panalty_26 + panalty_27 + panalty_28 + panalty_29...
         + panalty_30  + panalty_31 + panalty_32 + panalty_33 + panalty_34 + panalty_35 + panalty_36 + panalty_37 + panalty_38;      
end

%% Process Synthesis and Design Problems		

function o = Get_cec20_F8(x)
    %% Process synthesis problem
    [ps,D]=size(x);
    x(:,2) = round(x(:,2));
    f = 2*x(:,1) + x(:,2);
    g1 = 1.25 - x(:,1).^2 - x(:,2);
    g2 = x(:,1) + x(:,2) - 1.6;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    o = f + panalty_1 + panalty_2; 
end

function o = Get_cec20_F9(x)
    %% Process synthesis and design problem
    [ps,D]=size(x);
    x(:,3) = round(x(:,3));
    f = -x(:,3) + 2*x(:,1) + x(:,2);
    h1 = x(:,1) - 2*exp(-x(:,2));
    g1 = -x(:,1)+x(:,2)+x(:,3);
    H1 =abs(h1)-0.0001;
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,g1))^2;
    o = f + panalty_1 + panalty_2; 
end

function o = Get_cec20_F10(x)
    %% Process flow sheeting problem
    [ps,D]=size(x);
    x(:,3) = round(x(:,3));
    f = -0.7*x(:,3) + 5*(x(:,1)-0.5).^2 + 0.8;
    g1 = -exp(x(:,1) - 0.2) - x(:,2);
    g2 = x(:,2) + 1.1*x(:,3) + 1;
    g3 = x(:,1) - x(:,3) -0.2;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    o = f + panalty_1 + panalty_2 + panalty_3;
end

function o = Get_cec20_F11(x)
    %% Two-reactor Problem
    [ps,D]=size(x);
        x1 = x(:,1);
        x2 = x(:,2);
        v1 = x(:,3);
        v2 = x(:,4);
        y1 = round(x(:,5));
        y2 = round(x(:,6));
        x_ = x(:,7);
        
        z1 = 0.9*(1-exp(-0.5.*v1)).*x1;
        z2 = 0.8*(1-exp(-0.4*v2)).*x2;
        
        
        f = 7.5.*y1 + 5.5.*y2 + 7.*v1 + 6.*v2 + 5.*x_;
        
        h1 = y1 + y2 - 1;
        h2 = z1 + z2 - 10;
        h3 = x1 + x2 -x_;
        h4 = z1.*y1 + z2.*y2 - 10;
        g1 = v1 - 10*y1;
        g2 = v2 - 10*y2;
        g3 = x1 - 20*y1;
        g4 = x2 - 20*y2;
        panalty_1 = 10e100*(max(0,g1))^2;
        panalty_2 = 10e100*(max(0,g2))^2;
        panalty_3 = 10e100*(max(0,g3))^2;
        panalty_4 = 10e100*(max(0,g4))^2;
        H1 =abs(h1)-0.0001;
        H2 =abs(h2)-0.0001;
        H3 =abs(h3)-0.0001;
        H4 =abs(h4)-0.0001;
        panalty_5 = 10e100*(max(0,H1))^2;
        panalty_6 = 10e100*(max(0,H2))^2;
        panalty_7 = 10e100*(max(0,H3))^2;
        panalty_8 = 10e100*(max(0,H4))^2;
        o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8;
end

function o = Get_cec20_F12(x)
    %% Process synthesis problem
    [ps,D]=size(x);
        x1 = x(:,1);
        x2 = x(:,2);
        x3 = x(:,3);
        y1 = round(x(:,4));
        y2 = round(x(:,5));
        y3 = round(x(:,6));
        y4 = round(x(:,7));
        f = (y1-1).^2 + (y2-1).^2 + (y3-1).^2 - log(y4+1) + (x1-1).^22 + (x2-2).^2 + (x3-3).^2;
        g1 = x1 + x2 + x3 + y1 + y2 + y3 - 5;
        g2 = y3.^2 + x1.^2 + x2.^2 + x3.^2 - 5.5;
        g3 = x1 + y1 - 1.2;
        g4 = x2 + y2 - 1.8;
        g5 = x3 + y3 - 2.5;
        g6 = x1 + y4 - 1.2;
        g7 = y2.^2 + x2.^2 - 1.64;
        g8 = y3.^2 + x3.^2 - 4.25;
        g9 = y2.^2 + x3.^2 - 4.64;
        panalty_1 = 10e100*(max(0,g1))^2;
        panalty_2 = 10e100*(max(0,g2))^2;
        panalty_3 = 10e100*(max(0,g3))^2;
        panalty_4 = 10e100*(max(0,g4))^2;
        panalty_5 = 10e100*(max(0,g5))^2;
        panalty_6 = 10e100*(max(0,g6))^2;
        panalty_7 = 10e100*(max(0,g7))^2;
        panalty_8 = 10e100*(max(0,g8))^2;
        panalty_9 = 10e100*(max(0,g9))^2;
        o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9;
end

function o = Get_cec20_F13(x)
    %% Process design Problem
    [ps,D]=size(x);
        x1 = x(:,1);
        x2 = x(:,2);
        x3 = x(:,3);
        y1 = round(x(:,4));
        y2 = round(x(:,5));
        f = -5.357854*x1.^2 - 0.835689*y1.*x3 - 37.29329*y1 + 40792.141;
        a = [85.334407,0.0056858,0.0006262,0.0022053,80.51249, 0.0071317,....
                                                  0.0029955,0.0021813,9.300961,0.0047026,0.0012547,0.0019085];
        g1 = a(1) + a(2)*y2.*x3 + a(3)*y1.*x2 - a(4)*y1.*y1.*x3 - 92;
        g2 = a(5) + a(6)*y2.*x3 + a(7)*y1.*x2 + a(8)*x1.^2 -90 -20;
        g3 = a(9) + a(10)*y1.*x2 + a(11)*y1.*x1 + a(12)*x1.*x2 - 20 - 5;
        panalty_1 = 10e100*(max(0,g1))^2;
        panalty_2 = 10e100*(max(0,g2))^2;
        panalty_3 = 10e100*(max(0,g3))^2;
        o = f + panalty_1 + panalty_2 + panalty_3;
end


function o = Get_cec20_F14(x)
    %% Multi-product batch plant
    [ps,D]=size(x);
    %% constant
   S = [2,3,4;
        4,6,3];
   t = [8,20,8;
        16,4,4];
   H = 6000; alp = 250; beta = 0.6;
   Q1 = 40000; Q2 = 20000;
   %% decision Variable
   N1 = round(x(:,1)); N2 = round(x(:,2)); N3 = round(x(:,3));
   V1 = x(:,4); V2 = x(:,5); V3 = x(:,6);
   TL1 = x(:,7); TL2 = x(:,8);
   B1 = x(:,9); B2 = x(:,10);
   %% objective function
   f = alp.*(N1.*V1.^beta+N2.*V2.^beta+N3.*V3.^beta);
   %% constraints
   g1 = Q1.*TL1./B1+Q2.*TL2./B2-H;
   g2 = S(1,1).*B1+S(2,1).*B2-V1;
   g3 = S(1,2).*B1+S(2,2).*B2-V2;
   g4 = S(1,3).*B1+S(2,3).*B2-V3;
   g5 = t(1,1)-N1.*TL1;
   g6 = t(1,2)-N2.*TL1;
   g7 = t(1,3)-N3.*TL1;
   g8 = t(2,1)-N1.*TL2;
   g9 = t(2,2)-N2.*TL2;
   g10 = t(2,3)-N3.*TL2;
   panalty_1 = 10e100*(max(0,g1))^2;
   panalty_2 = 10e100*(max(0,g2))^2;
   panalty_3 = 10e100*(max(0,g3))^2;
   panalty_4 = 10e100*(max(0,g4))^2;
   panalty_5 = 10e100*(max(0,g5))^2;
   panalty_6 = 10e100*(max(0,g6))^2;
   panalty_7 = 10e100*(max(0,g7))^2;
   panalty_8 = 10e100*(max(0,g8))^2;
   panalty_9 = 10e100*(max(0,g9))^2;
   panalty_10 = 10e100*(max(0,g10))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9 + panalty_10;
end

%% Mechanical Engineering Problem		

function o = Get_cec20_F15(x)
    %% Weight Minimization of a Speed Reducer
    [ps,D]=size(x);
    f = 0.7854*x(:,1).*x(:,2).^2.*(3.3333.*x(:,3).^2+14.9334.*x(:,3)-43.0934)-1.508.*x(:,1).*(x(:,6).^2+x(:,7).^2).....
        +7.477.*(x(:,6).^3+x(:,7).^3)+0.7854.*(x(:,4).*x(:,6).^2+x(:,5).*x(:,7).^2);
  
    g1 = -x(:,1).*x(:,2).^2.*x(:,3)+27;
    g2 = -x(:,1).*x(:,2).^2.*x(:,3).^2+397.5;
    g3 = -x(:,2).*x(:,6).^4.*x(:,3).*x(:,4).^(-3)+1.93;
    g4 = -x(:,2).*x(:,7).^4.*x(:,3)./x(:,5).^3+1.93;
    g5 = 10.*x(:,6).^(-3).*sqrt(16.91.*10^6+(745.*x(:,4)./(x(:,2).*x(:,3))).^2)-1100;
    g6 = 10.*x(:,7).^(-3).*sqrt(157.5.*10^6+(745.*x(:,5)./(x(:,2).*x(:,3))).^2)-850;
    g7 = x(:,2).*x(:,3)-40;
    g8 = -x(:,1)./x(:,2)+5;
    g9 = x(:,1)./x(:,2)-12;
    g10 = 1.5.*x(:,6)-x(:,4)+1.9;
    g11 = 1.1.*x(:,7)-x(:,5)+1.9;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
    panalty_5 = 10e100*(max(0,g5))^2;
    panalty_6 = 10e100*(max(0,g6))^2;
    panalty_7 = 10e100*(max(0,g7))^2;
    panalty_8 = 10e100*(max(0,g8))^2;
    panalty_9 = 10e100*(max(0,g9))^2;
    panalty_10 = 10e100*(max(0,g10))^2;
    panalty_11 = 10e100*(max(0,g11))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9 + panalty_10 + panalty_11;
end

function o = Get_cec20_F16(x)
    %% Optimal Design of Industrial refrigeration System
    [ps,D]=size(x);
    f = 63098.88.*x(:,2).*x(:,4).*x(:,12)+5441.5.*x(:,2).^2.*x(:,12)+115055.5.*x(:,2).^1.664.*x(:,6).....
        +6172.27.*x(:,2).^2.*x(:,6)+63098.88.*x(:,1).*x(:,3).*x(:,11)+5441.5.*x(:,1).^2.*x(:,11).....
        +115055.5.*x(:,1).^1.664.*x(:,5)+6172.27.*x(:,1).^2.*x(:,5)+140.53.*x(:,1).*x(:,11)+281.29.*x(:,3).*x(:,11)....
        +70.26.*x(:,1).^2+281.29.*x(:,1).*x(:,3)+281.29.*x(:,3).^2+14437.*x(:,8).^1.8812.*x(:,12).^0.3424....
        .*x(:,10).*x(:,14).^(-1).*x(:,1).^2.*x(:,7).*x(:,9).^(-1)+20470.2.*x(:,7).^(2.893).*x(:,11).^0.316.*x(:,1).^2;
    g1 = 1.524.*x(:,7).^(-1)-1;
    g2 = 1.524.*x(:,8).^(-1)-1;
    g3 = 0.07789.*x(:,1)-2.*x(:,7).^(-1).*x(:,9)-1;
    g4 = 7.05305.*x(:,9).^(-1).*x(:,1).^2.*x(:,10).*x(:,8).^(-1).*x(:,2).^(-1).*x(:,14).^(-1)-1;
    g5 = 0.0833./x(:,13).*x(:,14)-1;
    g6 = 0.04771.*x(:,10).*x(:,8).^1.8812.*x(:,12).^0.3424-1;
    g7 = 0.0488.*x(:,9).*x(:,7).^1.893.*x(:,11).^0.316-1;
    g8 = 0.0099.*x(:,1)./x(:,3)-1;
    g9 = 0.0193.*x(:,2)./x(:,4)-1;
    g10 = 0.0298.*x(:,1)./x(:,5)-1;
    g11 = 47.136.*x(:,2).^0.333./x(:,10).*x(:,12)-1.333.*x(:,8).*x(:,13).^2.1195+62.08.*x(:,13).^2.1195.*x(:,8).^0.2./(x(:,12).*x(:,10))-1;
    g12 = 0.056.*x(:,2)./x(:,6)-1;
    g13 = 2./x(:,9)-1;
    g14 = 2./x(:,10)-1;
    g15 = x(:,12)./x(:,11)-1;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
    panalty_5 = 10e100*(max(0,g5))^2;
    panalty_6 = 10e100*(max(0,g6))^2;
    panalty_7 = 10e100*(max(0,g7))^2;
    panalty_8 = 10e100*(max(0,g8))^2;
    panalty_9 = 10e100*(max(0,g9))^2;
    panalty_10 = 10e100*(max(0,g10))^2;
    panalty_11 = 10e100*(max(0,g11))^2;
    panalty_12 = 10e100*(max(0,g12))^2;
    panalty_13 = 10e100*(max(0,g13))^2;
    panalty_14 = 10e100*(max(0,g14))^2;
    panalty_15 = 10e100*(max(0,g15))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9 + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15;
end

function o = Get_cec20_F17(x)
    %% Tension/compression  spring  design (case 1)
    [ps,D]=size(x);
    f = x(:,1).^2.*x(:,2).*(x(:,3)+2);
    g1 = 1-(x(:,2).^3.*x(:,3))./(71785.*x(:,1).^4);
    g2 = (4.*x(:,2).^2-x(:,1).*x(:,2))./(12566.*(x(:,2).*x(:,1).^3-x(:,1).^4))....
             + 1./(5108.*x(:,1).^2)-1;
    g3 = 1-140.45.*x(:,1)./(x(:,2).^2.*x(:,3));
    g4 = (x(:,1)+x(:,2))./1.5-1;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4;    
end

function o = Get_cec20_F18(x)
    %% update
    x(:,1) = 0.0625.*round(x(:,1));
    x(:,2) = 0.0625.*round(x(:,2));
    
    %% Pressure vessel design
    f = 0.6224.*x(:,1).*x(:,3).*x(:,4)+1.7781.*x(:,2).*x(:,3).^2....
        +3.1661.*x(:,1).^2.*x(:,4)+19.84.*x(:,1).^2.*x(:,3);
    g1 = -x(:,1)+0.0193.*x(:,3);
    g2 = -x(:,2)+0.00954.*x(:,3);
    g3 = -pi.*x(:,3).^2.*x(:,4)-4/3.*pi.*x(:,3).^3+1296000;
    g4 = x(:,4)-240;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4; 
end

function o = Get_cec20_F19(x)
    %% Welded beam design
    [ps,D]=size(x);
    f = 1.10471.*x(:,1).^2.*x(:,2)+0.04811.*x(:,3).*x(:,4).*(14+x(:,2));

    P = 6000; L = 14; delta_max = 0.25; E = 30*1e6; G = 12*1e6; 
    T_max = 13600; sigma_max = 30000;
    Pc = 4.013.*E.*sqrt(x(:,3).^2.*x(:,4).^6/30)./L.^2.*(1-x(:,3)./(2*L).*sqrt(E/(4.*G)));
    sigma = 6.*P.*L./(x(:,4).*x(:,3).^2);
    delta = 6.*P.*L.^3./(E.*x(:,3).^2.*x(:,4));
    J = 2.*(sqrt(2).*x(:,1).*x(:,2).*(x(:,2).^2/4+(x(:,1)+x(:,3)).^2./4));
    R = sqrt(x(:,2).^2/4+(x(:,1)+x(:,3)).^2/4);
    M = P.*(L+x(:,2)/2);
    ttt = M.*R./J;
    tt = P./(sqrt(2).*x(:,1).*x(:,2));
    t  = sqrt(tt.^2+2.*tt.*ttt.*x(:,2)./(2.*R)+ttt.^2);
    %% constraints
    g1 = t-T_max;
    g2 = sigma-sigma_max;
    g3 = x(:,1)-x(:,4);
    g4 = delta-delta_max;
    g5 = P-Pc;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
    panalty_5 = 10e100*(max(0,g5))^2;    
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5;     
end

function o = Get_cec20_F20(x)
    %% Three-bar truss design problem
    [ps,D]=size(x);
    f = (2.*sqrt(2).*x(:,1)+x(:,2))*100;
    g1 = (sqrt(2).*x(:,1)+x(:,2))./(sqrt(2).*x(:,1).^2+2.*x(:,1).*x(:,2))*2-2;
    g2 = x(:,2)./(sqrt(2).*x(:,1).^2+2.*x(:,1).*x(:,2))*2-2;
    g3 = 1./(sqrt(2).*x(:,2)+x(:,1))*2-2;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;   
    o = f + panalty_1 + panalty_2 + panalty_3;
end

function o = Get_cec20_F21(x)
    %% Multiple disk clutch brake design problem
    %% parameters
    [ps,D]=size(x);
   Mf = 3; Ms = 40; Iz = 55; n = 250; Tmax = 15; s = 1.5; delta = 0.5; 
   Vsrmax = 10; rho = 0.0000078; pmax = 1; mu = 0.6; Lmax = 30; delR = 20;
   Rsr = 2./3.*(x(:,2).^3-x(:,1).^3)./(x(:,2).^2.*x(:,1).^2);
   Vsr = pi.*Rsr.*n./30;
   A   = pi.*(x(:,2).^2-x(:,1).^2);
   Prz = x(:,4)./A;
   w   = pi.*n./30;
   Mh  = 2/3.*mu.*x(:,4).*x(:,5).*(x(:,2).^3-x(:,1).^3)./(x(:,2).^2-x(:,1).^2);
   T   = Iz.*w./(Mh+Mf);
   %%
   f = pi.*(x(:,2).^2-x(:,1).^2).*x(:,3).*(x(:,5)+1).*rho;
   g1 = -x(:,2)+x(:,1)+delR;
   g2 = (x(:,5)+1).*(x(:,3)+delta)-Lmax;
   g3 = Prz-pmax;
   g4 = Prz.*Vsr-pmax.*Vsrmax;
   g5 = Vsr-Vsrmax;
   g6 = T-Tmax;
   g7 = s.*Ms-Mh;
   g8 = -T;
   panalty_1 = 10e100*(max(0,g1))^2;
   panalty_2 = 10e100*(max(0,g2))^2;
   panalty_3 = 10e100*(max(0,g3))^2;
   panalty_4 = 10e100*(max(0,g4))^2;
   panalty_5 = 10e100*(max(0,g5))^2;
   panalty_6 = 10e100*(max(0,g6))^2;
   panalty_7 = 10e100*(max(0,g7))^2;
   panalty_8 = 10e100*(max(0,g8))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8; 
end

function o = Get_cec20_F22(x)
    %% Planetary gear train design optimization problem
    %% parameter Initialization
    [ps,D]=size(x);
   x = round(abs(x)); Pind = [3,4,5]; mind = [ 1.75, 2, 2.25, 2.5, 2.75, 3.0];
   N1 = x(:,1); N2 = x(:,2); N3 = x(:,3); N4 = x(:,4); N5 = x(:,5); N6 = x(:,6);
   p  = Pind(:,x(:,7))'; 
   m1 = mind(:,x(:,8))';
   m2 = mind(:,x(:,9))';
    %% objective function
    i1 = N6./N4; i01 = 3.11;
    i2 = N6.*(N1.*N3+N2.*N4)./(N1.*N3.*(N6-N4)); i02 = 1.84;
    iR = -(N2.*N6./(N1.*N3)); i0R = -3.11;
    f  = max([i1-i01,i2-i02,iR-i0R],[],2);
    %% constraints
    Dmax = 220; dlt22 = 0.5; dlt33 = 0.5; dlt55 = 0.5; dlt35 = 0.5; dlt34 = 0.5; dlt56 = 0.5;
    beta = acos(((N6-N3).^2+(N4+N5).^2-(N3+N5).^2)./(2.*(N6-N3).*(N4+N5)));
    g1 = m2.*(N6+2.5)-Dmax;
    g2 = m1.*(N1+N2)+m1.*(N2+2)-Dmax;
    g3 = m2.*(N4+N5)+m2.*(N5+2)-Dmax;
    g4 = abs(m1.*(N1+N2)-m2.*(N6-N3))-m1-m2;
    g5 = -((N1+N2).*sin(pi./p)-N2-2-dlt22);
    g6 = -((N6-N3).*sin(pi./p)-N3-2-dlt33);
    g7 = -((N4+N5).*sin(pi./p)-N5-2-dlt55);
    if beta == real(beta)
       g8 = (N3+N5+2+dlt35).^2-((N6-N3).^2+(N4+N5).^2-2.*(N6-N3).*(N4+N5).*cos(2.*pi./p-beta));
    else
       g8 = 1e6;
    end
    g9 = -(N6-2.*N3-N4-4-2.*dlt34);
    g10 = -(N6-N4-2.*N5-4-2.*dlt56);
    h1  = rem(N6-N4,p);
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
    panalty_5 = 10e100*(max(0,g5))^2;
    panalty_6 = 10e100*(max(0,g6))^2;
    panalty_7 = 10e100*(max(0,g7))^2;
    panalty_8 = 10e100*(max(0,g8))^2;
    panalty_9 = 10e100*(max(0,g9))^2;
    panalty_10 = 10e100*(max(0,g10))^2;
    H1 =abs(h1)-0.0001;
    panalty_11 = 10e100*(max(0,H1))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9 + panalty_10 + panalty_11;   
end

function o = Get_cec20_F23(x)
    %% Step-cone pulley problem
    %% parameter Initialization
    [ps,D]=size(x);
    d1 = x(:,1)*1e-3; d2 = x(:,2)*1e-3; d3 = x(:,3)*1e-3; d4 = x(:,4)*1e-3; w = x(:,5)*1e-3;
    N = 350; N1 = 750; N2 = 450; N3 = 250; N4 = 150;
    rho = 7200; a = 3; mu = 0.35; s = 1.75*1e6; t = 8*1e-3;
    %% objective function
    f = rho.*w.*pi./4.*(d1.^2.*(1+(N1./N).^2)+d2.^2.*(1+(N2./N).^2)+d3.^2.*(1+(N3./N).^2)+d4.^2.*(1+(N4./N).^2));
    %% constraint
    C1 = pi.*d1./2.*(1+N1./N)+(N1./N-1).^2.*d1.^2./(4.*a)+2.*a;
    C2 = pi.*d2./2.*(1+N2./N)+(N2./N-1).^2.*d2.^2./(4.*a)+2.*a;
    C3 = pi.*d3./2.*(1+N3./N)+(N3./N-1).^2.*d3.^2./(4.*a)+2.*a;
    C4 = pi.*d4./2.*(1+N4./N)+(N4./N-1).^2.*d4.^2./(4.*a)+2.*a;
    R1 = exp(mu.*(pi-2.*asin((N1./N-1).*d1./(2.*a))));
    R2 = exp(mu.*(pi-2.*asin((N2./N-1).*d2./(2.*a))));
    R3 = exp(mu.*(pi-2.*asin((N3./N-1).*d3./(2.*a))));
    R4 = exp(mu.*(pi-2.*asin((N4./N-1).*d4./(2.*a))));
    P1 = s.*t.*w.*(1-exp(-mu.*(pi-2.*asin((N1/N-1).*d1/(2.*a))))).*pi.*d1.*N1./60;
    P2 = s.*t.*w.*(1-exp(-mu.*(pi-2.*asin((N2/N-1).*d2/(2.*a))))).*pi.*d2.*N2./60;
    P3 = s.*t.*w.*(1-exp(-mu.*(pi-2.*asin((N3/N-1).*d3/(2.*a))))).*pi.*d3.*N3./60;
    P4 = s.*t.*w.*(1-exp(-mu.*(pi-2.*asin((N4/N-1).*d4/(2.*a))))).*pi.*d4.*N4./60;
    
    g1 = -R1+2;
    g2 = -R2+2;
    g3 = -R3+2;
    g4 = -R4+2;
    g5 = -P1+(0.75*745.6998);
    g6 = -P2+(0.75*745.6998);
    g7 = -P3+(0.75*745.6998);
    g8 = -P4+(0.75*745.6998);
    h1 = C1-C2;
    h2 = C1-C3;
    h3 = C1-C4;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
    panalty_5 = 10e100*(max(0,g5))^2;
    panalty_6 = 10e100*(max(0,g6))^2;
    panalty_7 = 10e100*(max(0,g7))^2;
    panalty_8 = 10e100*(max(0,g8))^2;
    H1 =abs(h1)-0.0001;
    panalty_9 = 10e100*(max(0,H1))^2;
    H2 =abs(h2)-0.0001;
    panalty_10 = 10e100*(max(0,H2))^2;
    H3 =abs(h3)-0.0001;
    panalty_11 = 10e100*(max(0,H3))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9 + panalty_10 + panalty_11;     
end

function o = Get_cec20_F24(x)
    %% Robot gripper problem
   [ps,D]=size(x);
   a = x(:,1); b = x(:,2); c = x(:,3); e = x(:,4); ff = x(:,5); l = x(:,6); delta = x(:,7);
   Ymin = 50; Ymax = 100; YG = 150; Zmax = 99.9999; P = 100;
   alpha_0 = acos((a.^2+l.^2+e.^2-b.^2)./(2.*a.*sqrt(l.^2+e.^2)))+atan(e./l);
   beta_0  = acos((b.^2+l.^2+e.^2-a.^2)./(2.*b.*sqrt(l.^2+e.^2)))-atan(e./l);
   alpha_m = acos((a.^2+(l-Zmax).^2+e.^2-b.^2)./(2.*a.*sqrt((l-Zmax).^2+e.^2)))+atan(e./(l-Zmax));
   beta_m  = acos((b.^2+(l-Zmax).^2+e.^2-a.^2)./(2.*b.*sqrt((l-Zmax).^2+e.^2)))-atan(e./(l-Zmax));
   %% objective function
   for i = 1:ps
       f(i,1) = -OBJ11(x(i,:),2)-OBJ11(x(i,:),1);
   end
   %% constraints
   Yxmin = 2.*(e+ff+c.*sin(beta_m+delta));
   Yxmax = 2.*(e+ff+c.*sin(beta_0+delta));
   g(:,1) = Yxmin-Ymin;
   g(:,2) = -Yxmin;
   g(:,3) = Ymax-Yxmax;
   g(:,4) = Yxmax-YG;
   g(:,5) = l.^2+e.^2-(a+b).^2;
   g(:,6) = b.^2-(a-e).^2-(l-Zmax).^2;
   g(:,7) = Zmax-l;

   tt     = imag(f)~=0;
   f(tt)  = 1e4;
   tt     = imag(g)~=0;
   g(tt)  = 1e4;
   panalty_1 = 10e100*(max(0,g(:,1)))^2;
   panalty_2 = 10e100*(max(0,g(:,2)))^2;
   panalty_3 = 10e100*(max(0,g(:,3)))^2;
   panalty_4 = 10e100*(max(0,g(:,4)))^2;
   panalty_5 = 10e100*(max(0,g(:,5)))^2;
   panalty_6 = 10e100*(max(0,g(:,6)))^2;
   panalty_7 = 10e100*(max(0,g(:,7)))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7;   
end

function o = Get_cec20_F25(x)
    %% Hydro-static thrust bearing design problem
    [ps,D]=size(x);
    R = x(:,1); Ro = x(:,2);  mu = x(:,3); Q = x(:,4);
    gamma = 0.0307; C = 0.5; n = -3.55; C1 = 10.04;
    Ws = 101000; Pmax = 1000; delTmax = 50; hmin = 0.001;
    gg = 386.4; N = 750;
    P    = (log10(log10(8.122*1e6.*mu+0.8))-C1)./n;
    delT = 2.*(10.^P-560);
    Ef   = 9336.*Q.*gamma.*C.*delT;
    h    = (2.*pi.*N./60).^2.*2.*pi.*mu./Ef.*(R.^4./4-Ro.^4./4)-1e-5;
    Po   = (6.*mu.*Q./(pi.*h.^3)).*log(R./Ro);
    W    = pi.*Po./2.*(R.^2-Ro.^2)./(log(R./Ro)-1e-5);
    %%  objective function
    f = (Q.*Po./0.7+Ef)./12;
    %%  constraints
    g1 = Ws-W;
    g2 = Po-Pmax;
    g3 = delT-delTmax;
    g4 = hmin-h;
    g5 = Ro-R;
    g6 = gamma./(gg.*Po).*(Q./(2.*pi.*R.*h))-0.001;
    g7 = W./(pi.*(R.^2-Ro.^2)+1e-5)-5000;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
    panalty_5 = 10e100*(max(0,g5))^2;
    panalty_6 = 10e100*(max(0,g6))^2;
    panalty_7 = 10e100*(max(0,g7))^2;
    o = f(1,1) + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7;  
end

function o = Get_cec20_F26(x)
    %% Four-stage gear box problem
    %% parameter initialized
    [ps,D]=size(x);
    x = round(x);
    Np1 = x(:,1); Ng1 = x(:,2); Np2 = x(:,3); Ng2 = x(:,4);
    Np3 = x(:,5); Ng3 = x(:,6); Np4 = x(:,7); Ng4 = x(:,8);
    Pvalue = [ 3.175, 5.715, 8.255, 12.7];
    b1 = Pvalue(x(:,9))'; b2 = Pvalue(x(:,10))'; b3 = Pvalue(x(:,11))';
    b4 = Pvalue(x(:,12))';
    XYvalue = [ 12.7,25.4,38.1,50.8,63.5,76.2,88.9,101.6,114.3];
    xp1 = XYvalue(x(:,13))'; xg1 = XYvalue(x(:,14))';
    xg2 = XYvalue(x(:,15))'; xg3 = XYvalue(x(:,16))';
    xg4 = XYvalue(x(:,17))'; yp1 = XYvalue(x(:,18))';
    yg1 = XYvalue(x(:,19))'; yg2 = XYvalue(x(:,20))';
    yg3 = XYvalue(x(:,21))'; yg4 = XYvalue(x(:,22))';
    %% value initilized
    c1 = sqrt((xg1-xp1).^2+(yg1-yp1).^2);
    c2 = sqrt((xg2-xp1).^2+(yg2-yp1).^2);
    c3 = sqrt((xg3-xp1).^2+(yg3-yp1).^2);
    c4 = sqrt((xg4-xp1).^2+(yg4-yp1).^2);
    CRmin = 1.4; dmin = 25.4; phi = 20.*pi./180; W = 55.9; JR = 0.2; Km = 1.6; Ko = 1.5; Lmax = 127; 
    sigma_H = 3290; sigma_N = 2090; w1 = 5000; wmin = 245; wmax = 255; Cp = 464;
    %% objective function
    f = pi./1000.*(b1.*c1.^2.*(Np1.^2+Ng1.^2)./(Np1+Ng1).^2+b2.*c2.^2.*(Np2.^2+Ng2.^2)./(Np2+Ng2).^2......
        +b3.*c3.^2.*(Np3.^2+Ng3.^2)./(Np3+Ng3).^2+b4.*c4.^2.*(Np4.^2+Ng4.^2)./(Np4+Ng4).^2);
    %% constraints
    g1 = (366000./(pi.*w1)+2.*c1.*Np1./(Np1+Ng1)).*((Np1+Ng1).^2./(4.*b1.*c1.^2.*Np1))-sigma_N.*JR./(0.0167.*W.*Ko.*Km);
    g2 = (366000.*Ng1./(pi.*w1.*Np1)+2.*c2.*Np2./(Np2+Ng2)).*((Np2+Ng2).^2./(4.*b2.*c2.^2.*Np2))-sigma_N.*JR./(0.0167.*W.*Ko.*Km);
    g3 = (366000.*Ng1.*Ng2./(pi.*w1.*Np1.*Np2)+2.*c3.*Np3./(Np3+Ng3)).*((Np3+Ng3).^2./(4.*b3.*c3.^2.*Np3))-sigma_N.*JR./(0.0167.*W.*Ko.*Km);
    g4 = (366000.*Ng1.*Ng2.*Ng3./(pi.*w1.*Np1.*Np2.*Np3)+2.*c4.*Np4./(Np4+Ng4)).*((Np4+Ng4).^2./(4.*b4.*c4.^2.*Np4))-sigma_N.*JR./(0.0167.*W.*Ko.*Km);
    g5 = (366000./(pi.*w1)+2.*c1.*Np1./(Np1+Ng1)).*((Np1+Ng1).^3./(4.*b1.*c1.^2.*Ng1.*Np1.^2))-(sigma_H./Cp).^2.*(sin(phi).*cos(phi))./(0.0334.*W.*Ko.*Km);
    g6 = (366000.*Ng1./(pi.*w1.*Np1)+2.*c2.*Np2./(Np2+Ng2)).*((Np2+Ng2).^3./(4.*b2.*c2.^2.*Ng2.*Np2.^2))-(sigma_H./Cp).^2.*(sin(phi).*cos(phi))./(0.0334.*W.*Ko.*Km);
    g7 = (366000.*Ng1.*Ng2./(pi.*w1.*Np1.*Np2)+2.*c3.*Np3./(Np3+Ng3)).*((Np3+Ng3).^3./(4.*b3.*c3.^2.*Ng3.*Np3.^2))-(sigma_H./Cp).^2.*(sin(phi).*cos(phi))./(0.0334.*W.*Ko.*Km);
    g8 = (366000.*Ng1.*Ng2.*Ng3./(pi.*w1.*Np1.*Np2.*Np3)+2.*c4.*Np4./(Np4+Ng4)).*((Np4+Ng4).^3./(4.*b4.*c4.^2.*Ng4.*Np4.^2))-(sigma_H./Cp).^2.*(sin(phi).*cos(phi))./(0.0334.*W.*Ko.*Km);
    g9 = CRmin.*pi.*cos(phi) - Np1.*sqrt(sin(phi).^2./4+1./Np1+(1./Np1).^2)-Ng1.*sqrt(sin(phi).^2./4+1./Ng1+(1./Ng1).^2)+sin(phi).*(Np1+Ng1)./2;
    g10 = CRmin.*pi.*cos(phi) - Np2.*sqrt(sin(phi).^2./4+1./Np2+(1./Np2).^2)-Ng2.*sqrt(sin(phi).^2./4+1./Ng2+(1./Ng2).^2)+sin(phi).*(Np2+Ng2)./2;
    g11 = CRmin.*pi.*cos(phi) - Np3.*sqrt(sin(phi).^2./4+1./Np3+(1./Np3).^2)-Ng3.*sqrt(sin(phi).^2./4+1./Ng3+(1./Ng3).^2)+sin(phi).*(Np3+Ng3)./2;
    g12 = CRmin.*pi.*cos(phi) - Np4.*sqrt(sin(phi).^2./4+1./Np4+(1./Np4).^2)-Ng4.*sqrt(sin(phi).^2./4+1./Ng4+(1./Ng4).^2)+sin(phi).*(Np4+Ng4)./2;
    g13 = dmin - 2.*c1.*Np1./(Np1+Ng1);
    g14 = dmin - 2.*c2.*Np2./(Np2+Ng2);
    g15 = dmin - 2.*c3.*Np3./(Np3+Ng3);
    g16 = dmin - 2.*c4.*Np4./(Np4+Ng4);
    g17 = dmin - 2.*c1.*Ng1./(Np1+Ng1);
    g18 = dmin - 2.*c2.*Ng2./(Np2+Ng2);
    g19 = dmin - 2.*c3.*Ng3./(Np3+Ng3);
    g20 = dmin - 2.*c4.*Ng4./(Np4+Ng4); 
    g21 = xp1 +((Np1+2).*c1./(Np1+Ng1))-Lmax;
    g22 = xg2 +((Np2+2).*c2./(Np2+Ng2))-Lmax;
    g23 = xg3 +((Np3+2).*c3./(Np3+Ng3))-Lmax;
    g24 = xg4 +((Np4+2).*c4./(Np4+Ng4))-Lmax;
    g25 = -xp1 +((Np1+2).*c1./(Np1+Ng1));
    g26 = -xg2 +((Np2+2).*c2./(Np2+Ng2));
    g27 = -xg3 +((Np3+2).*c3./(Np3+Ng3));
    g28 = -xg4 +((Np4+2).*c4./(Np4+Ng4));
    g29 = yp1 +((Np1+2).*c1./(Np1+Ng1))-Lmax;
    g30 = yg2 +((Np2+2).*c2./(Np2+Ng2))-Lmax;
    g31 = yg3 +((Np3+2).*c3./(Np3+Ng3))-Lmax;
    g32 = yg4 +((Np4+2).*c4./(Np4+Ng4))-Lmax;
    g33 = -yp1 +((Np1+2).*c1./(Np1+Ng1));
    g34 = -yg2 +((Np2+2).*c2./(Np2+Ng2));
    g35 = -yg3 +((Np3+2).*c3./(Np3+Ng3));
    g36 = -yg4 +((Np4+2).*c4./(Np4+Ng4));
    g37 = xg1 +((Ng1+2).*c1./(Np1+Ng1))-Lmax;
    g38 = xg2 +((Ng2+2).*c2./(Np2+Ng2))-Lmax;
    g39 = xg3 +((Ng3+2).*c3./(Np3+Ng3))-Lmax;
    g40 = xg4 +((Ng4+2).*c4./(Np4+Ng4))-Lmax;
    g41 = -xg1 +((Ng1+2).*c1./(Np1+Ng1));
    g42 = -xg2 +((Ng2+2).*c2./(Np2+Ng2));
    g43 = -xg3 +((Ng3+2).*c3./(Np3+Ng3));
    g44 = -xg4 +((Ng4+2).*c4./(Np4+Ng4));
    g45 = yg1 +((Ng1+2).*c1./(Np1+Ng1))-Lmax;
    g46 = yg2 +((Ng2+2).*c2./(Np2+Ng2))-Lmax;
    g47 = yg3 +((Ng3+2).*c3./(Np3+Ng3))-Lmax;
    g48 = yg4 +((Ng4+2).*c4./(Np4+Ng4))-Lmax;
    g49 = -yg1 +((Ng1+2).*c1./(Np1+Ng1));
    g50 = -yg2 +((Ng2+2).*c2./(Np2+Ng2));
    g51 = -yg3 +((Ng3+2).*c3./(Np3+Ng3));
    g52 = -yg4 +((Ng4+2).*c4./(Np4+Ng4));
    g53 = (0.945.*c1-Np1-Ng1).*(b1-5.715).*(b1-8.255).*(b1-12.70).*(-1);
    g54 = (0.945.*c2-Np2-Ng2).*(b2-5.715).*(b2-8.255).*(b2-12.70).*(-1);
    g55 = (0.945.*c3-Np3-Ng3).*(b3-5.715).*(b3-8.255).*(b3-12.70).*(-1);
    g56 = (0.945.*c4-Np4-Ng4).*(b4-5.715).*(b4-8.255).*(b4-12.70).*(-1);
    g57 = (0.646.*c1-Np1-Ng1).*(b1-3.175).*(b1-8.255).*(b1-12.70).*(+1);
    g58 = (0.646.*c2-Np2-Ng2).*(b2-3.175).*(b2-8.255).*(b2-12.70).*(+1);
    g59 = (0.646.*c3-Np3-Ng3).*(b3-3.175).*(b3-8.255).*(b3-12.70).*(+1);
    g60 = (0.646.*c4-Np4-Ng4).*(b4-3.175).*(b4-8.255).*(b4-12.70).*(+1);
    g61 = (0.504.*c1-Np1-Ng1).*(b1-3.175).*(b1-5.715).*(b1-12.70).*(-1);
    g62 = (0.504.*c2-Np2-Ng2).*(b2-3.175).*(b2-5.715).*(b2-12.70).*(-1);
    g63 = (0.504.*c3-Np3-Ng3).*(b3-3.175).*(b3-5.715).*(b3-12.70).*(-1);
    g64 = (0.504.*c4-Np4-Ng4).*(b4-3.175).*(b4-5.715).*(b4-12.70).*(-1);
    g65 = (0.0.*c1-Np1-Ng1).*(b1-3.175).*(b1-5.715).*(b1-8.255).*(+1);
    g66 = (0.0.*c2-Np2-Ng2).*(b2-3.175).*(b2-5.715).*(b2-8.255).*(+1);
    g67 = (0.0.*c3-Np3-Ng3).*(b3-3.175).*(b3-5.715).*(b3-8.255).*(+1);
    g68 = (0.0.*c4-Np4-Ng4).*(b4-3.175).*(b4-5.715).*(b4-8.255).*(+1);
    g69 = (-1.812.*c1+Np1+Ng1).*(b1-5.715).*(b1-8.255).*(b1-12.70).*(-1);
    g70 = (-1.812.*c2+Np2+Ng2).*(b2-5.715).*(b2-8.255).*(b2-12.70).*(-1);
    g71 = (-1.812.*c3+Np3+Ng3).*(b3-5.715).*(b3-8.255).*(b3-12.70).*(-1);
    g72 = (-1.812.*c4+Np4+Ng4).*(b4-5.715).*(b4-8.255).*(b4-12.70).*(-1);
    g73 = (-0.945.*c1+Np1+Ng1).*(b1-3.175).*(b1-8.255).*(b1-12.70).*(+1);
    g74 = (-0.945.*c2+Np2+Ng2).*(b2-3.175).*(b2-8.255).*(b2-12.70).*(+1);
    g75 = (-0.945.*c3+Np3+Ng3).*(b3-3.175).*(b3-8.255).*(b3-12.70).*(+1);
    g76 = (-0.945.*c4+Np4+Ng4).*(b4-3.175).*(b4-8.255).*(b4-12.70).*(+1);
    g77 = (-0.646.*c1+Np1+Ng1).*(b1-3.175).*(b1-5.715).*(b1-12.70).*(-1);
    g78 = (-0.646.*c2+Np2+Ng2).*(b2-3.175).*(b2-5.715).*(b2-12.70).*(-1);
    g79 = (-0.646.*c3+Np2+Ng3).*(b3-3.175).*(b3-5.715).*(b3-12.70).*(-1);
    g80 = (-0.646.*c4+Np3+Ng4).*(b4-3.175).*(b4-5.715).*(b4-12.70).*(-1);
    g81 = (-0.504.*c1+Np1+Ng1).*(b1-3.175).*(b1-5.715).*(b1-8.255).*(+1);
    g82 = (-0.504.*c2+Np2+Ng2).*(b2-3.175).*(b2-5.715).*(b2-8.255).*(+1);
    g83 = (-0.504.*c3+Np3+Ng3).*(b3-3.175).*(b3-5.715).*(b3-8.255).*(+1);
    g84 = (-0.504.*c4+Np4+Ng4).*(b4-3.175).*(b4-5.715).*(b4-8.255).*(+1);
    g85 = wmin -w1.*(Np1.*Np2.*Np3.*Np4)./(Ng1.*Ng2.*Ng3.*Ng4);
    g86 = -wmax +w1.*(Np1.*Np2.*Np3.*Np4)./(Ng1.*Ng2.*Ng3.*Ng4);
    g = [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, ...
     g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, ...
     g21, g22, g23, g24, g25, g26, g27, g28, g29, g30, ...
     g31, g32, g33, g34, g35, g36, g37, g38, g39, g40, ...
     g41, g42, g43, g44, g45, g46, g47, g48, g49, g50, ...
     g51, g52, g53, g54, g55, g56, g57, g58, g59, g60, ...
     g61, g62, g63, g64, g65, g66, g67, g68, g69, g70, ...
     g71, g72, g73, g74, g75, g76, g77, g78, g79, g80, ...
     g81, g82, g83, g84, g85, g86];
    g(g==inf) = 1e6;
    g(g==-inf) = 1e6;
    panalty_1 = 10e100*(max(0,g1))^2;
    panalty_2 = 10e100*(max(0,g2))^2;
    panalty_3 = 10e100*(max(0,g3))^2;
    panalty_4 = 10e100*(max(0,g4))^2;
    panalty_5 = 10e100*(max(0,g5))^2;
    panalty_6 = 10e100*(max(0,g6))^2;
    panalty_7 = 10e100*(max(0,g7))^2;
    panalty_8 = 10e100*(max(0,g8))^2;
    panalty_9 = 10e100*(max(0,g9))^2;
    panalty_10 = 10e100*(max(0,g10))^2;
    panalty_11 = 10e100*(max(0,g11))^2;
    panalty_12 = 10e100*(max(0,g12))^2;
    panalty_13 = 10e100*(max(0,g13))^2;
    panalty_14 = 10e100*(max(0,g14))^2;
    panalty_15 = 10e100*(max(0,g15))^2;
    panalty_16 = 10e100*(max(0,g16))^2;
    panalty_17 = 10e100*(max(0,g17))^2;
    panalty_18 = 10e100*(max(0,g18))^2;
    panalty_19 = 10e100*(max(0,g19))^2;
    panalty_20 = 10e100*(max(0,g20))^2;
    panalty_21 = 10e100*(max(0,g21))^2;
    panalty_22 = 10e100*(max(0,g22))^2;
    panalty_23 = 10e100*(max(0,g23))^2;
    panalty_24 = 10e100*(max(0,g24))^2;
    panalty_25 = 10e100*(max(0,g25))^2;
    panalty_26 = 10e100*(max(0,g26))^2;
    panalty_27 = 10e100*(max(0,g27))^2;
    panalty_28 = 10e100*(max(0,g28))^2;
    panalty_29 = 10e100*(max(0,g29))^2;
    panalty_30 = 10e100*(max(0,g30))^2;
    panalty_31 = 10e100*(max(0,g31))^2;
    panalty_32 = 10e100*(max(0,g32))^2;
    panalty_33 = 10e100*(max(0,g33))^2;
    panalty_34 = 10e100*(max(0,g34))^2;
    panalty_35 = 10e100*(max(0,g35))^2;
    panalty_36 = 10e100*(max(0,g36))^2;
    panalty_37 = 10e100*(max(0,g37))^2;
    panalty_38 = 10e100*(max(0,g38))^2;
    panalty_39 = 10e100*(max(0,g39))^2;
    panalty_40 = 10e100*(max(0,g40))^2;
    panalty_41 = 10e100*(max(0,g41))^2;
    panalty_42 = 10e100*(max(0,g42))^2;
    panalty_43 = 10e100*(max(0,g43))^2;
    panalty_44 = 10e100*(max(0,g44))^2;
    panalty_45 = 10e100*(max(0,g45))^2;
    panalty_46 = 10e100*(max(0,g46))^2;
    panalty_47 = 10e100*(max(0,g47))^2;
    panalty_48 = 10e100*(max(0,g48))^2;
    panalty_49 = 10e100*(max(0,g49))^2;
    panalty_50 = 10e100*(max(0,g50))^2;
    panalty_51 = 10e100*(max(0,g51))^2;
    panalty_52 = 10e100*(max(0,g52))^2;
    panalty_53 = 10e100*(max(0,g53))^2;
    panalty_54 = 10e100*(max(0,g54))^2;
    panalty_55 = 10e100*(max(0,g55))^2;
    panalty_56 = 10e100*(max(0,g56))^2;
    panalty_57 = 10e100*(max(0,g57))^2;
    panalty_58 = 10e100*(max(0,g58))^2;
    panalty_59 = 10e100*(max(0,g59))^2;
    panalty_60 = 10e100*(max(0,g60))^2;
    panalty_61 = 10e100*(max(0,g61))^2;
    panalty_62 = 10e100*(max(0,g62))^2;
    panalty_63 = 10e100*(max(0,g63))^2;
    panalty_64 = 10e100*(max(0,g64))^2;
    panalty_65 = 10e100*(max(0,g65))^2;
    panalty_66 = 10e100*(max(0,g66))^2;
    panalty_67 = 10e100*(max(0,g67))^2;
    panalty_68 = 10e100*(max(0,g68))^2;
    panalty_69 = 10e100*(max(0,g69))^2;
    panalty_70 = 10e100*(max(0,g70))^2;
    panalty_71 = 10e100*(max(0,g71))^2;
    panalty_72 = 10e100*(max(0,g72))^2;
    panalty_73 = 10e100*(max(0,g73))^2;
    panalty_74 = 10e100*(max(0,g74))^2;
    panalty_75 = 10e100*(max(0,g75))^2;
    panalty_76 = 10e100*(max(0,g76))^2;
    panalty_77 = 10e100*(max(0,g77))^2;
    panalty_78 = 10e100*(max(0,g78))^2;
    panalty_79 = 10e100*(max(0,g79))^2;
    panalty_80 = 10e100*(max(0,g80))^2;
    panalty_81 = 10e100*(max(0,g81))^2;
    panalty_82 = 10e100*(max(0,g82))^2;
    panalty_83 = 10e100*(max(0,g83))^2;
    panalty_84 = 10e100*(max(0,g84))^2;
    panalty_85 = 10e100*(max(0,g85))^2;
    panalty_86 = 10e100*(max(0,g86))^2;  
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15 + panalty_16 + panalty_17 + panalty_18 + panalty_19...
         + panalty_20 + panalty_21 + panalty_22 + panalty_23 + panalty_24 + panalty_25 + panalty_26 + panalty_27 + panalty_28 + panalty_29...
         + panalty_30  + panalty_31 + panalty_32 + panalty_33 + panalty_34 + panalty_35 + panalty_36 + panalty_37 + panalty_38 + panalty_39 + panalty_40...
         + panalty_41 + panalty_42 + panalty_43 + panalty_44 + panalty_45 + panalty_46 + panalty_47 + panalty_48 + panalty_49 + panalty_50...
         + panalty_51 + panalty_52 + panalty_53 + panalty_54 + panalty_55 + panalty_56 + panalty_57 + panalty_58 + panalty_59 + panalty_60...
         + panalty_61 + panalty_62 + panalty_63 + panalty_64 + panalty_65 + panalty_66 + panalty_67 + panalty_68 + panalty_69 + panalty_70...
         + panalty_71 + panalty_72 + panalty_73 + panalty_74 + panalty_75 + panalty_76 + panalty_77 + panalty_78 + panalty_79 + panalty_80...
         + panalty_81 + panalty_82 + panalty_83 + panalty_84 + panalty_85 + panalty_86;  
end

function o = Get_cec20_F27(x)
    %% 10-bar truss design
    [ps,D]=size(x);
    f = zeros(ps,1);
    % g = zeros(ps,3);
    % h = zeros(ps,1);
    
    for i = 1:ps
      %% objective function 
      f(i,1) = function_fitness(x(i,:));
      %% constraint
      [c,~] = ConsBar10(x(i,:));
      g(i,:)  = c;
    end
    panalty_1 = 10e100*(max(0,g(:,1)))^2;
    panalty_2 = 10e100*(max(0,g(:,2)))^2;
    panalty_3 = 10e100*(max(0,g(:,3)))^2;
    o = f(1,1) + panalty_1 + panalty_2 + panalty_3;    
end

function o = Get_cec20_F28(x)
    %% Rolling element bearing
    [ps,D]=size(x);
   Dm = x(:,1); Db = x(:,2); Z = round(x(:,3)); fi = x(:,4); fo = x(:,5);
   KDmin = x(:,6); KDmax = x(:,7); eps = x(:,8); e = x(:,9); chi = x(:,10);
   D = 160; d = 90; Bw = 30; T = D-d-2.*Db;
   phi_o = 2.*pi-2.*acos((((D-d).*0.5-0.75.*T).^2+(0.5.*D-0.25.*T-Db).^2-(0.5.*d+0.25.*T).^2)./(2.*(0.5.*(D-d)-0.75.*T).*(0.5.*D-0.25.*T-Db)));
   gamma = Db./Dm;
   fc    = 37.91.*(1+(1.04.*((1-gamma)./(1+gamma)).^1.72.*(fi.*(2.*fo-1)./(fo.*(2.*fi-1))).^0.41).^(10/3)).^(-0.3).....
          .*(gamma.^0.3.*(1-gamma).^1.39./(1+gamma).^(1./3)).*(2.*fi./(2.*fi-1)).^0.41;
   %% objective function
   ind    = find(Db > 25.4);
   f      = fc.*Z.^(2/3).*Db.^(1.8);
   f(ind) = 3.647.*fc(ind).*Z(ind).^(2/3).*Db(ind).^1.4;
   %% constraint
   g1 = Z-1-phi_o./(2.*asin(Db./Dm));
   g2 = KDmin.*(D-d)-2.*Db;
   g3 = 2.*Db-KDmax.*(D-d);
   g4 = chi.*Bw -Db;
   g5 = 0.5.*(D+d)-Dm;
   g6 = Dm-(0.5+e).*(D+d);
   g7 = eps.*Db-0.5.*(D-Dm-Db);
   g8 = 0.515 - fi;
   g9 = 0.515 - fo;
   panalty_1 = 10e100*(max(0,g1))^2;
   panalty_2 = 10e100*(max(0,g2))^2;
   panalty_3 = 10e100*(max(0,g3))^2;
   panalty_4 = 10e100*(max(0,g4))^2;
   panalty_5 = 10e100*(max(0,g5))^2;
   panalty_6 = 10e100*(max(0,g6))^2;
   panalty_7 = 10e100*(max(0,g7))^2;
   panalty_8 = 10e100*(max(0,g8))^2;
   panalty_9 = 10e100*(max(0,g9))^2;
   o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9;

end


function o = Get_cec20_F29(x)
    %% Gas Transmission Compressor Design (GTCD)
    [ps,D]=size(x);
    f = 8.61.*1e5.*x(:,1).^0.5.*x(:,2).*x(:,3).^(-2/3).*x(:,4).^(-1/2) +3.69.*1e4.*x(:,3)....
        + 7.72.*1e8.*x(:,1).^(-1).*x(:,2).^(0.219)-765.43.*1e6.*x(:,1).^(-1);
    g1 = x(:,4).*x(:,2).^(-2)+x(:,2).^(-2)-1;
    panalty_1 = 10e100*(max(0,g1))^2;
    o = f + panalty_1;    
end

function o = Get_cec20_F30(x)
    %% Tension/compression  spring  design (case 2)
    [ps,D]=size(x);
        x1 = round(x(:,1));
        x2 = x(:,2);
        d  = [0.009,0.0095,0.0104,0.0118,0.0128,0.0132,0.014,....
              0.015, 0.0162, 0.0173, 0.018, 0.020, 0.023, 0.025,...
              0.028, 0.032, 0.035, 0.041, 0.047, 0.054, 0.063,....
              0.072, 0.080, 0.092, 0.0105, 0.120, 0.135, 0.148,....
              0.162, 0.177, 0.192, 0.207, 0.225, 0.244, 0.263,....
              0.283, 0.307, 0.331, 0.362,0.394,0.4375,0.500];
        x3 = d(round(x(:,3))); x3 = x3(:);
       
        %% objective function
        f = (pi.^2.*x2.*x3.^2.*(x1+2))./4;
        %% constants
        cf = (4.*x2./x3-1)./(4.*x2./x3-4)+0.615.*x3./x2;
        K  = (11.5.*10.^6.*x3.^4)./(8.*x1.*x2.^3);
        lf = 1000./K + 1.05.*(x1+2).*x3;
        sigp = 300./K;
        g1 = (8000.*cf.*x2)./(pi.*x3.^3)-189000;
        g2 = lf-14;
        g3 = 0.2-x3;
        g4 = x2-3;
        g5 = 3-x2./x3;
        g6 = sigp - 6;
        g7 = sigp+700./K+1.05.*(x1+2).*x3-lf;
        g8 = 1.25-700./K;
        panalty_1 = 10e100*(max(0,g1))^2;
        panalty_2 = 10e100*(max(0,g2))^2;
        panalty_3 = 10e100*(max(0,g3))^2;
        panalty_4 = 10e100*(max(0,g4))^2;
        panalty_5 = 10e100*(max(0,g5))^2;
        panalty_6 = 10e100*(max(0,g6))^2;
        panalty_7 = 10e100*(max(0,g7))^2;
        panalty_8 = 10e100*(max(0,g8))^2;
        o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8;
end

function o = Get_cec20_F31(x)%//////////////////???????????????g and h///////////////////////////////
    %% Gear train  design Problem
    [ps,D]=size(x);
        x1 = x(:,1);
        x2 = x(:,2);
        x3 = x(:,3);
        x4 = x(:,4);
        f = (1/6.931-(x1.*x2)./(x3.*x4)).^2;
        g = zeros(ps,1);
        h = zeros(ps,1);
        o=f;
end

function o = Get_cec20_F32(x)
    %% Himmelblau's Function
    [ps,D]=size(x);
        x1 = x(:,1);
        x2 = x(:,2);
        x3 = x(:,3);
        x4 = x(:,4);
        x5 = x(:,5);
        %% objective function
        f = 5.3578547.*x3.^2 + 0.8356891.*x1.*x5 + 37.293239.*x1 - 40792.141;
        %% parameters
        G1 = 85.334407 + 0.0056858.*x2.*x5 + 0.0006262.*x1.*x4 - 0.0022053.*x3.*x5;
        G2 = 80.51249 + 0.0071317.*x2.*x5 + 0.0029955.*x1.*x2 + 0.0021813.*x3.^2;
        G3 = 9.300961 + 0.0047026.*x3.*x5 + 0.0012547.*x1.*x3 + 0.0019085.*x3.*x4;
        %% constraint
        g1 = G1-92;
        g2 = -G1;
        g3 = G2-110;
        g4 = -G2+90;
        g5 = G3-25;
        g6 = -G3+20;
        panalty_1 = 10e100*(max(0,g1))^2;
        panalty_2 = 10e100*(max(0,g2))^2;
        panalty_3 = 10e100*(max(0,g3))^2;
        panalty_4 = 10e100*(max(0,g4))^2;
        panalty_5 = 10e100*(max(0,g5))^2;
        panalty_6 = 10e100*(max(0,g6))^2;
        o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6;
end

function o = Get_cec20_F33(x)
    %% Topology Optimization
    [ps,D]=size(x);
     nely = 10;
     nelx = 3;
     penal = 3;
     for i = 1:ps
         X = [x(i,1:10);x(i,11:20);x(i,21:30)]';
         % FE-ANALYSIS
         [U]=FE(3,10,X,3);         
         % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
         [KE] = lk;
         c = 0.;
         for ely = 1:nely
             for elx = 1:nelx
                 n1 = (nely+1)*(elx-1)+ely; 
                 n2 = (nely+1)* elx   +ely;
                 Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
                 c = c + X(ely,elx)^penal*Ue'*KE*Ue;
                 dc(ely,elx) = -penal*X(ely,elx)^(penal-1)*Ue'*KE*Ue;
             end
         end
         % FILTERING OF SENSITIVITIES
         [dc]   = check(3,10,1.5,X,dc); 
         f(i,1) = c;
         g(i,:) = dc(1:end);
     end
     panalty=0;
     for ii=1:D
         panalty=panalty+10e100*(max(0,g(1,ii)))^2;
     end
     o = f + panalty;
end

%% Power System Problems		

function o = Get_cec20_F34(x)   
    global initial_flag
    persistent G B P Q
    %% Optimal Sizing of Single Phase Distributed Generation with reactive power support for Phase Balancing at Main Transformer/Grid
    [ps,D]=size(x);
    if initial_flag == 0
        disp('i entered');
        G = load('input data\FunctionPS1_G.txt');
        B = load('input data\FunctionPS1_B.txt');
        P = load('input data\FunctionPS1_P.txt');
        Q = load('input data\FunctionPS1_Q.txt');
        initial_flag = 1;
    end
    Y = G+1j*B;
    %% voltage initilization
    V = zeros(30,1);
    V(1)    = 1;
    V(2)    = exp(1j*4*pi/3);
    V(3)    = exp(1j*2*pi/3);
    Pdg     = zeros(30,1);
    Qdg     = zeros(30,1);
    Psp     = zeros(30,1);
    Qsp     = zeros(30,1);
    for i = 1:ps
        V(4:30)   = x(i,1:27)+1j*x(i,28:54);
        Psp(4:30) = x(i,55:81);
        Qsp(4:30) = x(i,82:108);
        Pdg([9,16,21,24,30]) = x(i,109:113);
        Qdg([9,16,21,24,30]) = x(i,114:118);
    %% current calculation
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delP = Psp-Pdg-P;
        delQ = Qsp-Qdg-Q;
        delIr = Ir-spIr;
        delIm = Im-spIm;
    %% objective calculation
       f(i,1) = abs(I(1)+I(2)+I(3))+abs(I(1)+exp(1j*4*pi/3)*I(2)+exp(1j*2*pi/3)*I(3));
       h(i,:) = [delIr(4:30)',delIm(4:30)',delP(4:30)',delQ(4:30)'];
    end
     panalty=0;
     for ii=1:108
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty;
end

function o = Get_cec20_F35(x)
    %% Optimal Sizing of Distributed Generation for Active Power Loss Minimization
    [ps,D]=size(x);
    global initial_flag
    persistent G B P Q
    if initial_flag == 0
        G = load('input data\FunctionPS2_G.txt');
        B = load('input data\FunctionPS2_B.txt');
        P = load('input data\FunctionPS2_P.txt');
        Q = load('input data\FunctionPS2_Q.txt');
        initial_flag = 1;
    end
    Y = G+1j*B;
    %% voltage initilization
    V = zeros(38,1);
    V(1)    = 1;
    Pdg     = zeros(38,1);
%     Qdg     = zeros(38,1);
    Psp     = zeros(38,1);
    Qsp     = zeros(38,1);
    for i = 1:ps
        V(2:38)   = x(i,1:37)+1j*x(i,38:74);
        Psp(2:38) = x(i,75:111);
        Qsp(2:38) = x(i,112:148);
        Pdg([34,35,36,37,38]) = x(i,149:153);
%         Qdg([34,35,36,37,38]) = x(i,154:158);
    %% current calculation
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delP = Psp-Pdg+P(:,1).*(abs(V)./P(:,5)).^P(:,6);
        delQ = Qsp+Q(:,1).*(abs(V)./Q(:,5)).^Q(:,6);
        delIr = Ir-spIr;
        delIm = Im-spIm;
    %% objective calculation
       f(i,1) = real(V(1).*conj(I(1)))+sum(Psp(2:38));
       h(i,:) = [delIr(2:38)',delIm(2:38)',delP(2:38)',delQ(2:38)'];
    end
     panalty=0;
     for ii=1:148
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty;
end

function o = Get_cec20_F36(x)
    %% Optimal Sizing of Distributed Generation (DG) and Capacitors for Reactive Power Loss Minimization
    [ps,D]=size(x);
    global initial_flag
    persistent G B P Q
    if initial_flag == 0
        G = load('input data\FunctionPS2_G.txt');
        B = load('input data\FunctionPS2_B.txt');
        P = load('input data\FunctionPS2_P.txt');
        Q = load('input data\FunctionPS2_Q.txt');
        initial_flag = 1;
    end
    Y = G+1j*B;
    %% voltage initilization
    V = zeros(38,1);
    V(1)    = 1;
    Pdg     = zeros(38,1);
    Qdg     = zeros(38,1);
    Psp     = zeros(38,1);
    Qsp     = zeros(38,1);
    for i = 1:ps
        V(2:38)   = x(i,1:37)+1j*x(i,38:74);
        Psp(2:38) = x(i,75:111);
        Qsp(2:38) = x(i,112:148);
        Pdg([34,35,36,37,38]) = x(i,149:153);
        Qdg([34,35,36,37,38]) = x(i,154:158);
    %% current calculation
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delP = Psp-Pdg+P(:,1).*(abs(V)./P(:,5)).^P(:,6);
        delQ = Qsp-Qdg+Q(:,1).*(abs(V)./Q(:,5)).^Q(:,6);
        delIr = Ir-spIr;
        delIm = Im-spIm;
    %% objective calculation
       f(i,1) = 0.5.*(real(V(1).*conj(I(1)))+sum(Psp(2:38)))+0.5.*(imag(V(1).*conj(I(1)))+sum(Qsp(2:38)));
       h(i,:) = [delIr(2:38)',delIm(2:38)',delP(2:38)',delQ(2:38)'];
    end
     panalty=0;
     for ii=1:148
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty;
end

function o = Get_cec20_F37(x)
    %% Optimal Power flow (Minimization of Active Power Loss)
    [ps,D]=size(x);
    global initial_flag
    persistent G B P Q
    if initial_flag == 0
        G = load('input data\FunctionPS11_G.txt');
        B = load('input data\FunctionPS11_B.txt');
        P = load('input data\FunctionPS11_P.txt');
        Q = load('input data\FunctionPS11_Q.txt');
        initial_flag = 1;
    end
    Y = G+1j*B;
    %% voltage initilization
    V = zeros(30,1);
    V(1)    = 1;
    Pg      = zeros(30,1);
    Qg      = zeros(30,1);
    Psp     = zeros(30,1);
    Qsp     = zeros(30,1);
    for i = 1:ps
        V(2:30)   = x(i,1:29)+ 1j*x(i,30:58);
        Psp(2:30) = x(i,59:87);
        Qsp(2:30) = x(i,88:116);
        Pg([2,13,22,23,27]) = x(i,117:121);
        Qg([2,13,22,23,27]) = x(i,122:126);
    %% current calculation
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delP = Psp-Pg+P;
        delQ = Qsp-Qg+Q;
        delIr = Ir-spIr;
        delIm = Im-spIm;
    %% objective calculation
       f(i,1) = real(V(1).*conj(I(1)))+sum(Psp(2:30));
       h(i,:) = [delIr(2:30)',delIm(2:30)',delP(2:30)',delQ(2:30)'];
    end
     panalty=0;
     for ii=1:116
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty;
end

function o = Get_cec20_F38(x)
    %% Optimal Power flow (Minimization of Fuel Cost)
    [ps,D]=size(x);
    global initial_flag
    persistent G B P Q
    if initial_flag == 0
        G = load('input data\FunctionPS11_G.txt');
        B = load('input data\FunctionPS11_B.txt');
        P = load('input data\FunctionPS11_P.txt');
        Q = load('input data\FunctionPS11_Q.txt');
        initial_flag = 1;
    end
    Y = G+1j*B;
    %% voltage initilization
    V = zeros(30,1);
    V(1)    = 1;
    Pg      = zeros(30,1);
    Qg      = zeros(30,1);
    Psp     = zeros(30,1);
    Qsp     = zeros(30,1);
    ng      = [1,2,13,22,23,27];
    a1      = [0,0,0,0,0,0];
    b1      = [ 2, 1.75,1,3.25,3,3];
    c1      = [ 0.02,0.0175,0.0625, 0.00834, 0.025, 0.0025];
    for i = 1:ps
        V(2:30)   = x(i,1:29)+ 1j*x(i,30:58);
        Psp(2:30) = x(i,59:87);
        Qsp(2:30) = x(i,88:116);
        Pg([2,13,22,23,27]) = x(i,117:121);
        Qg([2,13,22,23,27]) = x(i,122:126);
    %% current calculation
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delP = Psp-Pg+P;
        delQ = Qsp-Qg+Q;
        delIr = Ir-spIr;
        delIm = Im-spIm;
        Pg(1)= real(V(1).*conj(I(1)));
    %% objective calculation
       f(i,1) = sum(a1 + b1.*Pg(ng)' + c1.*(Pg(ng).^2)');
       h(i,:) = [delIr(2:30)',delIm(2:30)',delP(2:30)',delQ(2:30)'];
    end
     panalty=0;
     for ii=1:116
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty;
end

function o = Get_cec20_F39(x)
    %% Optimal Power flow (Minimization of Active Power Loss and Fuel Cost)
    [ps,D]=size(x);
    global initial_flag
    persistent G B P Q
    if initial_flag == 0
        G = load('input data\FunctionPS11_G.txt');
        B = load('input data\FunctionPS11_B.txt');
        P = load('input data\FunctionPS11_P.txt');
        Q = load('input data\FunctionPS11_Q.txt');
        initial_flag = 1;
    end
    Y = G+1j*B;
    %% voltage initilization
    V = zeros(30,1);
    V(1)    = 1;
    Pg      = zeros(30,1);
    Qg      = zeros(30,1);
    Psp     = zeros(30,1);
    Qsp     = zeros(30,1);
    ng      = [1,2,13,22,23,27];
    a1      = [0,0,0,0,0,0];
    b1      = [ 2, 1.75,1,3.25,3,3];
    c1      = [ 0.02,0.0175,0.0625, 0.00834, 0.025, 0.0025];
    for i = 1:ps
        V(2:30)   = x(i,1:29)+ 1j*x(i,30:58);
        Psp(2:30) = x(i,59:87);
        Qsp(2:30) = x(i,88:116);
        Pg([2,13,22,23,27]) = x(i,117:121);
        Qg([2,13,22,23,27]) = x(i,122:126);
    %% current calculation
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delP = Psp-Pg+P;
        delQ = Qsp-Qg+Q;
        delIr = Ir-spIr;
        delIm = Im-spIm;
        Pg(1)= real(V(1).*conj(I(1)));
    %% objective calculation
       f(i,1) = sum(a1 + b1.*Pg(ng)' + c1.*(Pg(ng).^2)')+0.75.*sum(Pg-P);
       h(i,:) = [delIr(2:30)',delIm(2:30)',delP(2:30)',delQ(2:30)'];
    end
     panalty=0;
     for ii=1:116
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty;
end

function o = Get_cec20_F40(x)
    %% Microgrid Power flow (Islanded case)
    [ps,D]=size(x);
    global initial_flag
    persistent P Q L
   if initial_flag == 0
        P = load('input data\FunctionPS2_P.txt');
        Q = load('input data\FunctionPS2_Q.txt');
        L = load('input data\FunctionPS14_linedata.txt');
        initial_flag = 1;
   end
    %% voltage initilization
    V = zeros(38,1);
    V(1)    = 1;
    Pc     = zeros(38,1);
    Qc     = zeros(38,1);
    for i = 1:ps
        V(2:38)   = x(i,1:37)+1j*x(i,38:74);
        Pc([34,35,36,37,38]) = 1./[5.102e-03;1.502e-03;4.506e-03;2.253e-03;2.253e-03];
        Qc([34,35,36,37,38]) = 1./[0.05;0.03;0.05;0.01;0.1];
        w         = x(i,75);
        V(1)      = x(i,76)+1e-5; 
    %% current calculation
        Y    = ybus(L,w);
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        Vr   = real(V);
        Vm   = imag(V);
        Psp  = Pc.*(1-w)-P(:,1).*(abs(V)./P(:,5)).^P(:,6);
        Qsp  = Qc.*(1-sqrt(Vr.^2+Vm.^2))-Q(:,1).*(abs(V)./Q(:,5)).^Q(:,6);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delIr = Ir-spIr;
        delIm = Im-spIm;
        delP2 = Psp-(Vr.*Ir+Vm.*Im);
        delQ2 = Qsp-(Vm.*Ir-Vr.*Im);
    %% objective calculation
       f(i,1) = sum(delP2(1:38).^2)+sum(delQ2(1:38).^2) ;
       h(i,:) = [delIr(1:38)',delIm(1:38)'];
    end
     panalty=0;
     for ii=1:76
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty;
end

function o = Get_cec20_F41(x)
    %% Microgrid Power flow (Grid-connected case)
    [ps,D]=size(x);
    global initial_flag
    persistent G B P Q
   if initial_flag == 0
        G = load('input data\FunctionPS2_G.txt');
        B = load('input data\FunctionPS2_B.txt');
        P = load('input data\FunctionPS2_P.txt');
        Q = load('input data\FunctionPS2_Q.txt');
        initial_flag = 1;
    end
    Y = G+1j*B;
    %% voltage initilization
    V       = zeros(38,1);
    V(1)    = 1;
    Pdg     = zeros(38,1);
    Qdg     = zeros(38,1);
    Pdg([34,35,36,37,38]) = 0.2;
    Qdg([34,35,36,37,38]) = 0.18;
    for i = 1:ps
        V(2:38)   = x(i,1:37)+1j*x(i,38:74);
    %% current calculation
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        Vr   = real(V);
        Vm   = imag(V);
        Psp = Pdg-P(:,1).*(abs(V)./P(:,5)).^P(:,6);
        Qsp = Qdg-Q(:,1).*(abs(V)./Q(:,5)).^Q(:,6);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delIr = Ir-spIr;
        delIm = Im-spIm;
        delP = Psp-(Vr.*Ir+Vm.*Im);
        delQ = Qsp-(Vm.*Ir-Vr.*Im);
    %% objective calculation
       f(i,1) = sum(delP(2:38).^2)+sum(delQ(2:38).^2) ;
       h(i,:) = [delIr(2:38)',delIm(2:38)'];
    end
     panalty=0;
     for ii=1:74
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty; 
end

function o = Get_cec20_F42(x)
    %% Optimal Setting of Droop Controller for Minimization of Active Power Loss in Islanded Microgrids
    [ps,D]=size(x);
    global initial_flag
    persistent P Q L
   if initial_flag == 0
        P = load('input data\FunctionPS2_P.txt');
        Q = load('input data\FunctionPS2_Q.txt');
        L = load('input data\FunctionPS14_linedata.txt');
        initial_flag = 1;
   end
    %% voltage initilization
    V = zeros(38,1);
    V(1)    = 1;
    Pc     = zeros(38,1);
    Qc     = zeros(38,1);
    for i = 1:ps
        V(2:38)   = x(i,1:37)+1j*x(i,38:74);
        w         = x(i,75);
        V(1)      = x(i,76)+1e-5;
        Pc([34,35,36,37,38]) = x(i,77:81);
        Qc([34,35,36,37,38]) = x(i,82:86);
    %% current calculation
        Y    = ybus(L,w);
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        Vr   = real(V);
        Vm   = imag(V);
        Psp  = Pc.*(1-w)-P(:,1).*(abs(V)./P(:,5)).^P(:,6);
        Qsp  = Qc.*(1-sqrt(Vr.^2+Vm.^2))-Q(:,1).*(abs(V)./Q(:,5)).^Q(:,6);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delIr = Ir-spIr;
        delIm = Im-spIm;
    %% objective calculation
       f(i,1) = sum(Psp) ;
       h(i,:) = [delIr(1:38)',delIm(1:38)'];
    end
     panalty=0;
     for ii=1:76
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty; 
end

function o = Get_cec20_F43(x)
    %% Optimal Setting of Droop Controller for Minimization of  Reactive Power Loss in Islanded Microgrids
    [ps,D]=size(x);
    global initial_flag
    persistent P Q L
    if initial_flag == 0
        P = load('input data\FunctionPS2_P.txt');
        Q = load('input data\FunctionPS2_Q.txt');
        L = load('input data\FunctionPS14_linedata.txt');
        initial_flag = 1;
   end
    %% voltage initilization
    V = zeros(38,1);
    V(1)    = 1;
    Pc     = zeros(38,1);
    Qc     = zeros(38,1);
    for i = 1:ps
        V(2:38)   = x(i,1:37)+1j*x(i,38:74);
        w         = x(i,75);
        V(1)      = x(i,76)+1e-5;
        Pc([34,35,36,37,38]) = x(i,77:81);
        Qc([34,35,36,37,38]) = x(i,82:86);
    %% current calculation
        Y    = ybus(L,w);
        I    = Y*V;
        Ir   = real(I);
        Im   = imag(I);
        Vr   = real(V);
        Vm   = imag(V);
        Psp  = Pc.*(1-w)-P(:,1).*(abs(V)./P(:,5)).^P(:,6);
        Qsp  = Qc.*(1-sqrt(Vr.^2+Vm.^2))-Q(:,1).*(abs(V)./Q(:,5)).^Q(:,6);
        spI  = conj((Psp+1j*Qsp)./V);
        spIr = real(spI);
        spIm = imag(spI);
        delIr = Ir-spIr;
        delIm = Im-spIm;
    %% objective calculation
       f(i,1) = 0.5*(sum(Qsp)+sum(Psp));
       h(i,:) = [delIr(1:38)',delIm(1:38)'];
    end
     panalty=0;
     for ii=1:76
         H1 =abs(h(1,ii))-0.0001;
         panalty2 = 10e100*(max(0,H1))^2;
         panalty=panalty+panalty2;
     end
     o = f + panalty; 
   
end

function o = Get_cec20_F44(x)
    %% Wind Farm Layout Problem
    [ps,D]=size(x);
    interval      = 15;                        
    interval_num  = fix(360 / interval);     
    cut_in_speed  = 3.5;                  
    rated_speed   = 14;                    
    cut_out_speed = 25;                   
    R             = 40;                               
    H             = 80;                               
    CT            = 0.8;                             
    a             = 1 - sqrt(1 - CT);                 
    kappa         = 0.01;                         
    minDistance   = 5 * R;                  
    N             = 15;                               
    X             = 2000;                             
    Y             = 2000;                             
    k(1 : interval_num) = 2;
    c = [7 5 5 5 5 4 5 6 7 7 8 9.5 10 8.5 8.5 6.5 4.6 2.6 8 5 6.4 5.2 4.5 3.9];
    fre = [0.0003	0.0072	0.0237	0.0242	0.0222	0.0301	0.0397	0.0268	0.0626 ...	
          0.0801	0.1025	0.1445	0.1909	0.1162	0.0793	0.0082	0.0041	0.0008 ...	
          0.0010	0.0005	0.0013	0.0031	0.0085	0.0222];
    %% Objective Function
    for i = 1:ps
        f(i,1) = - Fitness(interval_num, interval, fre, N,x(i,:), ...,
                   a, kappa, R, k, c, cut_in_speed, rated_speed, cut_out_speed, 'origin');
    end
    %% Constraint Violation
    for i = 1:(0.5*D)
        XX(:,i) = x(:,2*i-1);
        YY(:,i) = x(:,2*i);
    end
    k = 1;
    for i = 1:(0.5*D)
        for j = (i+1):(0.5*D)-1
           g(:,k) = 5*R - sqrt((XX(:,i)-XX(:,j)).^2+(YY(:,i)-YY(:,j)).^2);
           k = k+1;
        end
    end
     panalty=0;
     for ii=1:91
         panalty=panalty+10e100*(max(0,g(:,ii)))^2;
     end
     o = f + panalty; 
end

%% Power Electronic Problems		

function o = Get_cec20_F45(x)
    %% SOPWM for 3-level Invereters
    [ps,D]=size(x);
    m = 0.32;
    s = (-ones(1,25)).^(2:26);
    k = [5,7,11,13,17,19,23,25,29,31,35,37,41,43,47,49,53,55,59,61,65,67,71,73,77,79,83,85,91,95,97];
    for i = 1:ps
        su = 0;
        for j = 1:31
            su2 = 0;
            for l = 1:D
                su2 = su2 + s(l).*cos(k(j).*x(i,l)*pi/180);
            end
            su = su + su2.^2./k(j).^4;
        end
        f(i,1) = (su).^0.5./(sum(1./k.^4)).^0.5;
    end
    g = zeros(ps,D-1);
    for i = 1:D-1
        g(:,i) = x(:,i)-x(:,i+1)+1e-6;
    end
     panalty=0;
     for ii=1:24
         panalty=panalty+10e100*(max(0,g(:,ii)))^2;
     end
    h = sum(s.*cos(x*pi/180),2)-m;
    H1 =abs(h)-0.0001;
    panalty2 = 10e100*(max(0,H1))^2;    
    o = f + panalty + panalty2; 
end

function o = Get_cec20_F46(x)
    %% SOPWM for 5-level Inverters
    [ps,D]=size(x);
    m = 0.32;
    s = [1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,1,1,-1];
    k = [5,7,11,13,17,19,23,25,29,31,35,37,41,43,47,49,53,55,59,61,65,67,71,73,77,79,83,85,91,95,97];
    for i = 1:ps
        su = 0;
        for j = 1:31
            su2 = 0;
            for l = 1:D
                su2 = su2 + s(l).*cos(k(j).*x(i,l)*pi/180);
            end
            su = su + su2.^2./k(j).^4;
        end
        f(i,1) = 0.5.*(su).^0.5./(sum(1./k.^4)).^0.5;
    end
    g = zeros(ps,D-1);
    for i = 1:D-1
        g(:,i) = x(:,i)-x(:,i+1)+1e-6;
    end
     panalty=0;
     for ii=1:24
         panalty=panalty+10e100*(max(0,g(:,ii)))^2;
     end    
    h = sum(s.*cos(x*pi/180),2)-2*m;
    H1 =abs(h)-0.0001;
    panalty2 = 10e100*(max(0,H1))^2;    
    o = f + panalty + panalty2;     
end

function o = Get_cec20_F47(x)
    %% SOPWM for 7-level Inverters
    [ps,D]=size(x);
    m = 0.36;
    s = [1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,1];
    k = [5,7,11,13,17,19,23,25,29,31,35,37,41,43,47,49,53,55,59,61,65,67,71,73,77,79,83,85,91,95,97];
    for i = 1:ps
        su = 0;
        for j = 1:31
            su2 = 0;
            for l = 1:D
                su2 = su2 + s(l).*cos(k(j).*x(i,l)*pi/180);
            end
            su = su + su2.^2./k(j).^4;
        end
        f(i,1) = 1/3.*(su).^0.5./(sum(1./k.^4)).^0.5;
    end
    g = zeros(ps,D-1);
    for i = 1:D-1
        g(:,i) = x(:,i)-x(:,i+1)+1e-6;
    end
     panalty=0;
     for ii=1:24
         panalty=panalty+10e100*(max(0,g(:,ii)))^2;
     end       
    h = sum(s.*cos(x*pi/180),2)-3*m;   
    H1 =abs(h)-0.0001;
    panalty2 = 10e100*(max(0,H1))^2;    
    o = f + panalty + panalty2;      
end

function o = Get_cec20_F48(x)
    %% SOPWM for 9-level Inverters
    [ps,D]=size(x);
    m = 0.32;
    s = [1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,1,1,1,1,-1,1];
    k = [5,7,11,13,17,19,23,25,29,31,35,37,41,43,47,49,53,55,59,61,65,67,71,73,77,79,83,85,91,95,97];
    for i = 1:ps
        su = 0;
        for j = 1:31
            su2 = 0;
            for l = 1:D
                su2 = su2 + s(l).*cos(k(j).*x(i,l)*pi/180);
            end
            su = su + su2.^2./k(j).^4;
        end
        f(i,1) = 1/4.*(su).^0.5./(sum(1./k.^4)).^0.5;
    end
    g = zeros(ps,D-1);
    for i = 1:D-1
        g(:,i) = x(:,i)-x(:,i+1)+1e-6;
    end
     panalty=0;
     for ii=1:29
         panalty=panalty+10e100*(max(0,g(:,ii)))^2;
     end         
    h = sum(s.*cos(x*pi/180),2)-4*m; 
    H1 =abs(h)-0.0001;
    panalty2 = 10e100*(max(0,H1))^2;    
    o = f + panalty + panalty2;     
end

function o = Get_cec20_F49(x)
    %% SOPWM for 11-level Inverters
    [ps,D]=size(x);
    m = 0.3333;
    s = [1,-1,1,1,1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,1,1,-1,-1,1,-1,-1];
    k = [5,7,11,13,17,19,23,25,29,31,35,37,41,43,47,49,53,55,59,61,65,67,71,73,77,79,83,85,91,95,97];
    for i = 1:ps
        su = 0;
        for j = 1:31
            su2 = 0;
            for l = 1:D
                su2 = su2 + s(l).*cos(k(j).*x(i,l)*pi/180);
            end
            su = su + su2.^2./k(j).^4;
        end
        f(i,1) = 1/5.*(su).^0.5./(sum(1./k.^4)).^0.5;
    end
    g = zeros(ps,D-1);
    for i = 1:D-1
        g(:,i) = x(:,i)-x(:,i+1)+1e-6;
    end
     panalty=0;
     for ii=1:29
         panalty=panalty+10e100*(max(0,g(:,ii)))^2;
     end       
    h = sum(s.*cos(x*pi/180),2)-5*m; 
    H1 =abs(h)-0.0001;
    panalty2 = 10e100*(max(0,H1))^2;    
    o = f + panalty + panalty2;     
end

function o = Get_cec20_F50(x)
    %% SOPWM for 13-level Inverters
    [ps,D]=size(x);
    m = 0.32;
    s = [1,1,1,-1,1,-1,1,-1,1,1,1,1,-1,-1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1,-1,1,-1,1];
    k = [5,7,11,13,17,19,23,25,29,31,35,37,41,43,47,49,53,55,59,61,65,67,71,73,77,79,83,85,91,95,97];
    for i = 1:ps
        su = 0;
        for j = 1:31
            su2 = 0;
            for l = 1:D
                su2 = su2 + s(l).*cos(k(j).*x(i,l)*pi/180);
            end
            su = su + su2.^2./k(j).^4;
        end
        f(i,1) = 1/6.*(su).^0.5./(sum(1./k.^4)).^0.5;
    end
    g = zeros(ps,D-1);
    for i = 1:D-1
        g(:,i) = x(:,i)-x(:,i+1)+1e-6;
    end
     panalty=0;
     for ii=1:29
         panalty=panalty+10e100*(max(0,g(:,ii)))^2;
     end         
    h = sum(s.*cos(x*pi/180),2)-6*m; 
    H1 =abs(h)-0.0001;
    panalty2 = 10e100*(max(0,H1))^2;    
    o = f + panalty + panalty2;      
end

%% Livestock Feed Ration Optimization 		

function o = Get_cec20_F51(x)
    %% Livestock Feed Ration Optimization (case 1)
    [ps,D]=size(x);
    global initial_flag
    persistent P
    if initial_flag == 0
        P = load('input data\FunctionRM_feed.txt');
        initial_flag = 1;     
    end
    %% objective function
    f = sum(x.*repmat(P(1,:),ps,1),2);
    %% constraints
    g1 = -sum(x.*repmat(P(5,:),ps,1),2)+1.090;
    g2 = sum(x.*repmat(P(5,:),ps,1),2)-2.170;
    g3 = -sum(x.*repmat(P(4,:),ps,1),2)+4.870;
    g4 = sum(x.*repmat(P(4,:),ps,1),2)-5.200;
    g5 = -sum(x.*repmat(P(6,:),ps,1),2)+0.043;
    g6 = sum(x.*repmat(P(6,:),ps,1),2)-0.086; 
    g7 = -sum(x.*repmat(P(7,:),ps,1),2)+0.023; 
    g8 = sum(x.*repmat(P(7,:),ps,1),2)-0.046;
    g9 = -sum(x(:,1:17),2)./sum(x,2)+0.295; 
    g10= sum(x(:,1:17),2)./sum(x,2)-0.36; 
    g11= -sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2)+0.3; 
    g12= sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2)-0.4712; 
    g13= -sum(x(:,34:59),2)+9.2; 
    g14= sum(x(:,34:59),2)-11.5; 
    h1 = sum(x.*repmat(P(2,:),ps,1),2)-6.9; 
      panalty_1 = 10e100*(max(0,g1))^2;
      panalty_2 = 10e100*(max(0,g2))^2;
      panalty_3 = 10e100*(max(0,g3))^2;
      panalty_4 = 10e100*(max(0,g4))^2;
      panalty_5 = 10e100*(max(0,g5))^2;
      panalty_6 = 10e100*(max(0,g6))^2;
      panalty_7 = 10e100*(max(0,g7))^2;
      panalty_8 = 10e100*(max(0,g8))^2;
      panalty_9 = 10e100*(max(0,g9))^2;
      panalty_10 = 10e100*(max(0,g10))^2;
      panalty_11 = 10e100*(max(0,g11))^2;
      panalty_12 = 10e100*(max(0,g12))^2;
      panalty_13 = 10e100*(max(0,g13))^2;
      panalty_14 = 10e100*(max(0,g14))^2;
      H1 =abs(h1)-0.0001;
      panalty_15 = 10e100*(max(0,H1))^2; 
      o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15; 
end

function o = Get_cec20_F52(x)
    %% Livestock Feed Ration Optimization (case 2)
    [ps,D]=size(x);
    global initial_flag
    persistent P
    if initial_flag == 0
        P = load('input data\FunctionRM_feed.txt');
        initial_flag = 1;     
    end
    %% objective function
    f = sum(x.*repmat(P(1,:),ps,1),2); 
    %% constraint
    g1 = -sum(x.*repmat(P(5,:),ps,1),2)+1.280; 
    g2 = sum(x.*repmat(P(5,:),ps,1),2)-2.560 ; 
    g3 = -sum(x.*repmat(P(4,:),ps,1),2)+7.300; 
    g4 = sum(x.*repmat(P(4,:),ps,1),2)-7.810; 
    g5 = -sum(x.*repmat(P(6,:),ps,1),2)+0.005;
    g6 = sum(x.*repmat(P(6,:),ps,1),2)-0.094; 
    g7 = -sum(x.*repmat(P(7,:),ps,1),2)+0.031;
    g8 = sum(x.*repmat(P(7,:),ps,1),2)-0.062;
    g9 = -sum(x(:,1:17),2)./sum(x,2)+0.2;
    g10= sum(x(:,1:17),2)./sum(x,2)-0.24; 
    g11= -sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2)+0.3;
    g12= sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2)-0.4; 
    g13= -sum(x(:,34:59),2)+9.8; 
    g14 = sum(x(:,34:59),2)-16.4; 
    h1 = sum(x.*repmat(P(2,:),ps,1),2)-9.8; 
      panalty_1 = 10e100*(max(0,g1))^2;
      panalty_2 = 10e100*(max(0,g2))^2;
      panalty_3 = 10e100*(max(0,g3))^2;
      panalty_4 = 10e100*(max(0,g4))^2;
      panalty_5 = 10e100*(max(0,g5))^2;
      panalty_6 = 10e100*(max(0,g6))^2;
      panalty_7 = 10e100*(max(0,g7))^2;
      panalty_8 = 10e100*(max(0,g8))^2;
      panalty_9 = 10e100*(max(0,g9))^2;
      panalty_10 = 10e100*(max(0,g10))^2;
      panalty_11 = 10e100*(max(0,g11))^2;
      panalty_12 = 10e100*(max(0,g12))^2;
      panalty_13 = 10e100*(max(0,g13))^2;
      panalty_14 = 10e100*(max(0,g14))^2;
      H1 =abs(h1)-0.0001;
      panalty_15 = 10e100*(max(0,H1))^2; 
      o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15;     
end


function o = Get_cec20_F53(x)
    %% Livestock Feed Ration Optimization (case 3)
    [ps,D]=size(x);
    global initial_flag
    persistent P
    if initial_flag == 0
        P = load('input data\FunctionRM_feed.txt');
        initial_flag = 1;     
    end
    %% objective function
    f = sum(x.*repmat(P(1,:),ps,1),2); 
    %% constraint
    g1 = -sum(x.*repmat(P(5,:),ps,1),2)+1.170;
    g2 = sum(x.*repmat(P(5,:),ps,1),2)-2.340;
    g3 = -sum(x.*repmat(P(4,:),ps,1),2)+6.940; 
    g4 = sum(x.*repmat(P(4,:),ps,1),2)-7.430; 
    g5 = -sum(x.*repmat(P(6,:),ps,1),2)+0.038; 
    g6 = sum(x.*repmat(P(6,:),ps,1),2)-0.076;
    g7 = -sum(x.*repmat(P(7,:),ps,1),2)+0.034; 
    g8 = sum(x.*repmat(P(7,:),ps,1),2)-0.068;
    g9 = -sum(x(:,1:17),2)./sum(x,2)+0.085;
    g10= sum(x(:,1:17),2)./sum(x,2)-0.111;
    g11= -sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2)+0.25;
    g12= sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2)-0.4;
    g13= -sum(x(:,34:59),2)+11.6;
    g14= sum(x(:,34:59),2)-14.5; 
    h1 = sum(x.*repmat(P(2,:),ps,1),2)-8.7; 
      panalty_1 = 10e100*(max(0,g1))^2;
      panalty_2 = 10e100*(max(0,g2))^2;
      panalty_3 = 10e100*(max(0,g3))^2;
      panalty_4 = 10e100*(max(0,g4))^2;
      panalty_5 = 10e100*(max(0,g5))^2;
      panalty_6 = 10e100*(max(0,g6))^2;
      panalty_7 = 10e100*(max(0,g7))^2;
      panalty_8 = 10e100*(max(0,g8))^2;
      panalty_9 = 10e100*(max(0,g9))^2;
      panalty_10 = 10e100*(max(0,g10))^2;
      panalty_11 = 10e100*(max(0,g11))^2;
      panalty_12 = 10e100*(max(0,g12))^2;
      panalty_13 = 10e100*(max(0,g13))^2;
      panalty_14 = 10e100*(max(0,g14))^2;
      H1 =abs(h1)-0.0001;
      panalty_15 = 10e100*(max(0,H1))^2; 
      o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15;  
end

function o = Get_cec20_F54(x)
    %% Livestock Feed Ration Optimization (case 4)
    [ps,D]=size(x);
    global initial_flag
    persistent P    
    if initial_flag == 0
        P = load('input data\FunctionRM_feed.txt');
        initial_flag = 1;     
    end
    %% objective function
    f = sum(x.*repmat(P(1,:),ps,1),2); 
    %% constraints
    g1 = -sum(x.*repmat(P(5,:),ps,1),2) + 0.56 ; 
    g2 = sum(x.*repmat(P(5,:),ps,1),2) - 1.12 ; 
    g3 = -sum(x.*repmat(P(4,:),ps,1),2) + 3.23 ; 
    g4 = sum(x.*repmat(P(4,:),ps,1),2) - 3.46 ; 
    g5 = -sum(x.*repmat(P(6,:),ps,1),2) + 0.018 ;
    g6 = sum(x.*repmat(P(6,:),ps,1),2) - 0.036 ;
    g7 = -sum(x.*repmat(P(7,:),ps,1),2) + 0.0116 ;
    g8 = sum(x.*repmat(P(7,:),ps,1),2) - 0.040 ; 
    g9 = -sum(x(:,1:17),2)./sum(x,2) + 0.25 ; 
    g10 = sum(x(:,1:17),2)./sum(x,2) - 0.9 ; 
    g11 = -sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2) + 0.3; 
    g12 = sum(x.*repmat(P(3,:),ps,1),2)./sum(x,2) - 0.4384; 
    g13 = -sum(x(:,34:59),2) + 7.470; 
    g14 = sum(x(:,34:59),2) - 9.340; 
    h1 = sum(x.*repmat(P(2,:),ps,1),2)-5.6; 
      panalty_1 = 10e100*(max(0,g1))^2;
      panalty_2 = 10e100*(max(0,g2))^2;
      panalty_3 = 10e100*(max(0,g3))^2;
      panalty_4 = 10e100*(max(0,g4))^2;
      panalty_5 = 10e100*(max(0,g5))^2;
      panalty_6 = 10e100*(max(0,g6))^2;
      panalty_7 = 10e100*(max(0,g7))^2;
      panalty_8 = 10e100*(max(0,g8))^2;
      panalty_9 = 10e100*(max(0,g9))^2;
      panalty_10 = 10e100*(max(0,g10))^2;
      panalty_11 = 10e100*(max(0,g11))^2;
      panalty_12 = 10e100*(max(0,g12))^2;
      panalty_13 = 10e100*(max(0,g13))^2;
      panalty_14 = 10e100*(max(0,g14))^2;
      H1 =abs(h1)-0.0001;
      panalty_15 = 10e100*(max(0,H1))^2; 
      o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6 + panalty_7 + panalty_8 + panalty_9...
         + panalty_10 + panalty_11 + panalty_12 + panalty_13 + panalty_14 + panalty_15;  
end

function o = Get_cec20_F55(x)
    %% Livestock Feed Ration Optimization (case 5)
    [ps,D]=size(x);
    global initial_flag
    persistent P       
    if initial_flag == 0
        P = load('input data\FunctionRM_dairy.txt');
        initial_flag = 1;     
    end
    %% objective function
    f = sum(x.*repmat(P(1,:),ps,1),2); 
    %% constraints
    h1   =  sum(x.*repmat(P(12,:),ps,1),2)-25.67;
    h2   =  sum(x.*repmat(P(2,:),ps,1),2)-0.0218;
    h3   =  sum(x.*repmat(P(3,:),ps,1),2)-0.062;
    h4   =  sum(x.*repmat(P(13,:),ps,1),2)-0.034; 
    h5   =  sum(x.*repmat(P(14,:),ps,1),2)-0.021;
    h6   =  sum(repmat(sum(P(2:11,:),1),ps,1).*x,2)-0.999; 
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    H5 =abs(h5)-0.0001;
    H6 =abs(h6)-0.0001;
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;
    panalty_5 = 10e100*(max(0,H5))^2;
    panalty_6 = 10e100*(max(0,H6))^2;
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6;
end

function o = Get_cec20_F56(x)
    %% Livestock Feed Ration Optimization (case 6)
    [ps,D]=size(x);
    global initial_flag
    persistent P  
    if initial_flag == 0
        P = load('input data\FunctionRM_dairy.txt');
        initial_flag = 1;     
    end
    %% objective function
    f = sum(x.*repmat(P(1,:),ps,1),2); 
    %% constraints
    h1   =  sum(x.*repmat(P(12,:),ps,1),2)-65.24;
    h2   =  sum(x.*repmat(P(2,:),ps,1),2)-0.066; 
    h3   =  sum(x.*repmat(P(3,:),ps,1),2)-0.159; 
    h4   =  sum(x.*repmat(P(13,:),ps,1),2)-0.103;
    h5   =  sum(x.*repmat(P(14,:),ps,1),2)-0.052; 
    h6   =  sum(repmat(sum(P(2:11,:),1),ps,1).*x,2)-2.644; 
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    H5 =abs(h5)-0.0001;
    H6 =abs(h6)-0.0001;
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;
    panalty_5 = 10e100*(max(0,H5))^2;
    panalty_6 = 10e100*(max(0,H6))^2;
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6;
end

function o = Get_cec20_F57(x)
    %% Livestock Feed Ration Optimization (case 7)
    [ps,D]=size(x);
    global initial_flag
    persistent P  
    if initial_flag == 0
        P = load('input data\FunctionRM_dairy.txt');
        initial_flag = 1;     
    end
    %% objective function
    f = sum(x.*repmat(P(1,:),ps,1),2); 
    %% constraints
    h1   =  sum(x.*repmat(P(12,:),ps,1),2)-30.05; 
    h2   =  sum(x.*repmat(P(2,:),ps,1),2)-0.0259;
    h3   =  sum(x.*repmat(P(3,:),ps,1),2)-0.077;
    h4   =  sum(x.*repmat(P(13,:),ps,1),2)-0.096; 
    h5   =  sum(x.*repmat(P(14,:),ps,1),2)-0.025;
    h6   =  sum(repmat(sum(P(2:11,:),1),ps,1).*x,2)-1.214; 
    H1 =abs(h1)-0.0001;
    H2 =abs(h2)-0.0001;
    H3 =abs(h3)-0.0001;
    H4 =abs(h4)-0.0001;
    H5 =abs(h5)-0.0001;
    H6 =abs(h6)-0.0001;
    panalty_1 = 10e100*(max(0,H1))^2;
    panalty_2 = 10e100*(max(0,H2))^2;
    panalty_3 = 10e100*(max(0,H3))^2;
    panalty_4 = 10e100*(max(0,H4))^2;
    panalty_5 = 10e100*(max(0,H5))^2;
    panalty_6 = 10e100*(max(0,H6))^2;
    o = f + panalty_1 + panalty_2 + panalty_3 + panalty_4 + panalty_5 + panalty_6;
end


% g=g';
% h=h';
% 
% end

% Program to for Admittance And Impedance Bus Formation....

function Y = ybus(linedata,f)  % Returns Y
linedata(:,4) = linedata(:,4).*f;
% linedata(:,3:4) = linedata(:,3:4).*10000/127^2;
% linedata(:,3:4) = linedata(:,3:4);


fb = linedata(:,1);             % From bus number...
tb = linedata(:,2);             % To bus number...
r = linedata(:,3);              % Resistance, R...
x = linedata(:,4);              % Reactance, X...
b = linedata(:,5);              % Ground Admittance, B/2...
a = linedata(:,6);              % Tap setting value..
z = r + i*x;                    % z matrix...
y = 1./z;                       % To get inverse of each element...
b = i*b;                        % Make B imaginary...

nb = max(max(fb),max(tb));      % No. of buses...
nl = length(fb);                % No. of branches...
Y = zeros(nb,nb);               % Initialise YBus...

 % Formation of the Off Diagonal Elements...
 for k = 1:nl
     Y(fb(k),tb(k)) = Y(fb(k),tb(k)) - y(k)/a(k);
     Y(tb(k),fb(k)) = Y(fb(k),tb(k));
 end

 % Formation of Diagonal Elements....
 for m = 1:nb
     for n = 1:nl
         if fb(n) == m
             Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             Y(m,m) = Y(m,m) + y(n) + b(n);
         end
     end
 end
end
function ff = OBJ11(x,n)
a = x(1); b = x(2); c = x(3); e = x(4); f = x(5); l = x(6); 
 Zmax = 99.9999; P = 100;
if n == 1
     fhd = @(z) P.*b.*sin(acos((a.^2+(l-z).^2+e.^2-b.^2)./(2.*a.*sqrt((l-z).^2+e.^2)))+acos((b.^2+(l-z).^2+e.^2-a.^2)./(2.*b.*sqrt((l-z).^2+e.^2))))./....
       (2.*c.*cos(acos((a.^2+(l-z).^2+e.^2-b.^2)./(2.*a.*sqrt((l-z).^2+e.^2)))+atan(e./(l-z))));
else
    fhd = @(z) -(P.*b.*sin(acos((a.^2+(l-z).^2+e.^2-b.^2)./(2.*a.*sqrt((l-z).^2+e.^2)))+acos((b.^2+(l-z).^2+e.^2-a.^2)./(2.*b.*sqrt((l-z).^2+e.^2))))./....
       (2.*c.*cos(acos((a.^2+(l-z).^2+e.^2-b.^2)./(2.*a.*sqrt((l-z).^2+e.^2)))+atan(e./(l-z)))));
end
options = optimset('Display','off');
 [~,ff]= fminbnd(fhd,0,Zmax,options); 
end

function [Weight] = function_fitness(section)

E   = 6.98*1e10;      % Young's elastic modulus (N/m^2)
A   = section;        % area of bar (m^2)
rho = 2770;           % density of material (kg/m^3)
%--------------------------------------------------------------------------
%           1         2       3       4       5     6                     
gcoord = [18.288,  18.288,  9.144,  9.144,      0,  0 
           9.144,       0,  9.144,      0,  9.144,  0];
%          1  2  3  4  5  6  7  8  9  10
element = [3, 1, 4, 2, 3, 1, 4, 3, 2, 1
           5, 3, 6, 4, 4, 2, 5, 6, 3, 4];
%--------------------------------------------------------------------------
% calculate Weight matrix
Weight = 0;
for i=1:length(element)
    nd = element(:,i);
    x  = gcoord(1,nd); y = gcoord(2,nd);
    % compute long of each bar
    le = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    Weight =  Weight + rho*le*A(i);
end
end
function [c,ceq] = ConsBar10(x)
type = '2D';
E    = 6.98*1e10;      % Young's elastic modulus (N/m^2)
A    = x;
rho  = 2770;           % density of material (kg/m^3)
%--------------------------------------------------------------------------
%           1        2        3       4       5     6          
gcoord  = [18.288,  18.288,  9.144,  9.144,      0,  0 
           9.144,       0,  9.144,      0,  9.144,  0];
%          1  2  3  4  5  6  7  8  9  10
element = [3, 1, 4, 2, 3, 1, 4, 3, 2, 1
           5, 3, 6, 4, 4, 2, 5, 6, 3, 4];
nel     = length(element);    % total element
nnode   = length(gcoord);     % total node
ndof    = 2;                  % number of degree of freedom of one node
sdof    = nnode*ndof;         % total dgree of freedom of system
% plotModel( type,gcoord,element );
% calculate stiffness matrix 
[ K,M ] = Cal_K_and_M( type,gcoord,element,A,rho,E );
% add non-structural mass
addedMass = 454; %kg
for idof = 1:sdof
    M(idof,idof) = M(idof,idof) + addedMass;
end
% apply boundary
bcdof   = [(5:6)*2-1, (5:6)*2];     % boundary condition displacement
% Giai phuong trinh tim tri rieng va vector rieng
[omega_2,~]=eigens(K,M,bcdof); 
f=sqrt(omega_2)/2/pi;
% f(1:5)
c1 = 7/f(1) -1;
c2 = 15/f(2)-1;
c3 = 20/f(3)-1;
c = [c1,c2,c3];
ceq = [];
end

function [ K,M ] = Cal_K_and_M( type,gcoord,element,A,rho,E )
% calculate K and M
nel     = length(element);    % total element
nnode   = length(gcoord);     % total node
switch type
    case '3D'
        ndof    = 3;                  % number of degree of freedom of one node
        sdof    = nnode*ndof;         % total dgree of freedom of system
        K       = zeros(sdof,sdof);
        M       = zeros(sdof,sdof);
        for iel=1:nel
            nd = element(:,iel);
            x  = gcoord(1,nd); y = gcoord(2,nd); z = gcoord(3,nd);
            % compute long of each bar
            le = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2 + (z(2)-z(1))^2);
            % compute direction cosin
            l_ij = (x(2)-x(1))/le;      % Eq.8.19
            m_ij = (y(2)-y(1))/le;      % Eq.8.19
            n_ij = (z(2)-z(1))/le;      % Eq.8.19
            % compute transform matrix
            Te = [l_ij m_ij  n_ij   0       0     0;
                0    0      0   l_ij   m_ij   n_ij];
            % compute stiffness matrix of element
            ke = A(iel)*E/le*[1 -1; -1  1];
            ke = Te'*ke*Te;
            me = rho*le*A(iel)*[2 0 0 1 0 0
                0 2 0 0 1 0;
                0 0 2 0 0 1;
                1 0 0 2 0 0;
                0 1 0 0 2 0;
                0 0 1 0 0 2]/6;
            % find index assemble
            index   = [3*nd(1)-2 3*nd(1)-1 3*nd(1)  3*nd(2)-2 3*nd(2)-1  3*nd(2)];
            % assemble ke in K
            K(index,index) = K(index,index) + ke;
            M(index,index) = M(index,index) + me;
        end

    case '2D'
        ndof    = 2;                  % number of degree of freedom of one node
        sdof    = nnode*ndof;         % total dgree of freedom of system
        K       = zeros(sdof,sdof);
        M       = zeros(sdof,sdof);
        for iel=1:nel
            nd = element(:,iel);
            x  = gcoord(1,nd); y = gcoord(2,nd);
            % compute long of each bar
            le = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
            % compute direction cosin
            l_ij = (x(2)-x(1))/le;
            m_ij = (y(2)-y(1))/le;
            % compute transform matrix
            Te = [l_ij m_ij   0      0 ;
                0    0   l_ij   m_ij];

            % compute stiffness matrix of element
            ke = A(iel)*E/le*[1 -1;
                -1  1];
            ke = Te'*ke*Te;
            me = rho*le*A(iel)*[2 0 1 0;
                0 2 0 1
                1 0 2 0
                0 1 0 2]/6; % lumped mass matrix
            % find index assemble
            index   = [2*nd(1)-1 2*nd(1)  2*nd(2)-1  2*nd(2)];
            % assemble ke in K
            K(index,index) = K(index,index) + ke;
            % assemble me in M
            M(index,index) = M(index,index) + me;
        end
end
end

function [L,X]=eigens(K,M,b)
  [nd,nd]=size(K);
  fdof=[1:nd]';
%
  if nargin==3
    pdof=b(:);
    fdof(pdof)=[]; 
    if nargout==2
      [X1,D]=eig(K(fdof,fdof),M(fdof,fdof));
      [nfdof,nfdof]=size(X1);
      for j=1:nfdof;
        mnorm=sqrt(X1(:,j)'*M(fdof,fdof)*X1(:,j));
        X1(:,j)=X1(:,j)/mnorm;
      end
      d=diag(D);
      [L,i]=sort(d);
      X2=X1(:,i);
      X=zeros(nd,nfdof);
      X(fdof,:)=X2;
    else
      d=eig(K(fdof,fdof),M(fdof,fdof));
      L=sort(d);
    end
  else
    if nargout==2
      [X1,D]=eig(K,M);
      for j=1:nd;
        mnorm=sqrt(X1(:,j)'*M*X1(:,j));
        X1(:,j)=X1(:,j)/mnorm;
      end
      d=diag(D);
      [L,i]=sort(d);
      X=X1(:,i);
    else
      d=eig(K,M);
      L=sort(d);
    end
  end
end


function all_power = Fitness(interval_num, interval, fre, N, coordinate, ...,
            a, kappa, R, k, c, cut_in_speed, rated_speed, cut_out_speed, evaluate_method)
all_power = 0;                 
for i = 1 : interval_num
   interval_dir = (i - 0.5) * interval;
   [power_eva] = eva_power(i, interval_dir, N, coordinate, ...,
            a, kappa, R,k(i), c(i), cut_in_speed, rated_speed, cut_out_speed, evaluate_method);
    all_power = all_power + fre(i) * sum(power_eva);
end
end

function power_eva = eva_power(interval_dir_num, interval_dir, N, coordinate, ...,
           a, kappa, R, k, c, cut_in_speed, rated_speed, cut_out_speed, evaluate_method)

if(strcmp(evaluate_method, 'caching'))
    [vel_def] = eva_func_deficit_caching(interval_dir_num ,N, coordinate, interval_dir, a, kappa, R);
else
    [vel_def] = eva_func_deficit(interval_dir_num, N, coordinate, interval_dir, a, kappa, R);
end
interval_c(1 : N) = 0;
for i = 1 : N
   interval_c(i) = c * (1 - vel_def(i)); 
end
n_ws = (rated_speed - cut_in_speed) / 0.3;
power_eva(1 : N) = 0;
for i = 1 : N
    for j = 1 : n_ws
        v_j_1 = cut_in_speed + (j - 1) * 0.3;
        v_j = cut_in_speed + j * 0.3;
        power_eva(i) = power_eva(i) + 1500 * exp((v_j_1 + v_j) / 2 - 7.5) ./ (5 + exp((v_j_1 + v_j) / 2 - 7.5)) * ...,
            (exp(-(v_j_1 / interval_c(i))^k) - exp(-(v_j / interval_c(i))^k));
    end
    power_eva(i) = power_eva(i) + 1500 * (exp(-(rated_speed / interval_c(i))^k) - exp(-(cut_out_speed / interval_c(i))^k));
end
end


function[vel_def] = eva_func_deficit(interval_dir_num, N, coordinate, theta, a, kappa, R)


global thetaVeldefijMatrix;

vel_def(1 : N) = 0;

for i = 1 : N
    vel_def_i = 0;
    for j = 1 : N   
        [affected, dij] = downstream_wind_turbine_is_affected(coordinate, j, i, theta, kappa, R);
        if(affected)  
            def = a / (1 + kappa * dij / R)^2;
%             def = restrict(def, 1);
            thetaVeldefijMatrix(i, j, interval_dir_num) = def;
            vel_def_i = vel_def_i + def^2;  
        else
            thetaVeldefijMatrix(i, j, interval_dir_num) = 0;
        end  
    end
%     vel_def_i = restrict(vel_def_i, 1);
    vel_def(i) = sqrt(vel_def_i);
end
end

function[vel_def] = eva_func_deficit_caching(interval_dir_num, N, coordinate, theta, a, kappa, R)

global thetaVeldefijMatrix;
global turbineMoved;

vel_def(1 : N) = 0;
movedTurbine = 1;
for i = 1 : N
    if(turbineMoved(i) == 1)
        movedTurbine = i;
    end
end

for i = 1 : N

    vel_def_i = 0;

    if(i ~= movedTurbine)
        [affected, dij] = downstream_wind_turbine_is_affected(coordinate, movedTurbine, i, theta, kappa, R);
        if(affected)  
            def = a / (1 + kappa * dij / R)^2;
            def = restrict(def, 1);
        else      
            def = 0;
        end 
        vel_def_i = sum((thetaVeldefijMatrix(i, :, interval_dir_num)).^2) - (thetaVeldefijMatrix(i, movedTurbine, interval_dir_num))^2 + def^2;
        thetaVeldefijMatrix(i, movedTurbine, interval_dir_num) = def;
    else
        for j = 1 : N   
            [affected, dij] = downstream_wind_turbine_is_affected(coordinate, j, i, theta, kappa, R);
            if(affected)  
                def = a / (1 + kappa * dij / R)^2;
                def = restrict(def, 1);
            else
                def = 0;      
            end
            vel_def_i = vel_def_i + def^2; 
            thetaVeldefijMatrix(i,j,interval_dir_num) = def;
        end
    end
    vel_def_i = restrict(vel_def_i, 1);
    vel_def(i) = sqrt(vel_def_i);
end
end

function[affected, dij] = downstream_wind_turbine_is_affected(coordinate, upstream_wind_turbine, ...,
    downstream_wind_turbine, theta, kappa, R)

    affected = 0;
    Tijx = (coordinate(2 * downstream_wind_turbine - 1) - coordinate(2 * upstream_wind_turbine - 1));
    Tijy = (coordinate(2 * downstream_wind_turbine) - coordinate(2 * upstream_wind_turbine));
    dij = cosd(theta) * Tijx + sind(theta) * Tijy;
    lij = sqrt((Tijx^2 + Tijy^2) - (dij)^2);
    l = dij * kappa + R;
    if((upstream_wind_turbine ~= downstream_wind_turbine) && (l > lij-R) && (dij > 0))
        affected = 1;
    end
end



%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%F(2,1) = -10000;
F(2*(nely+1)*(nelx+1),1)=-10000; 
%fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
fixeddofs   = [1:2*(nely+1)];
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
end
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 206000000.; 
nu = 0.3;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end

