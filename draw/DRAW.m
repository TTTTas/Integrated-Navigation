clear all;
close all;
[Sow_KF, ...
    X_KF, Y_KF, Z_KF, ...
    B_KF, L_KF, H_KF, ...
    E_KF, N_KF, U_KF, ...
    clock_KF, ...
    Vx_KF, Vy_KF, Vz_KF, ...
    VC_KF,...
    GS_KF, BS_KF] = importKF("Static-KF.kf");
[Sow_LS, ...
    X_LS, Y_LS, Z_LS, ...
    B_LS, L_LS, H_LS, ...
    E_LS, N_LS, U_LS, ...
    clock_LS,...
    ~, ~, ...
    Vx_LS, Vy_LS, Vz_LS, ...
    VC_LS, ...
    ~, ~, ~,...
    GS_LS,BS_LS] = importLS("Static-LS.pos");
%% 
img = figure(WindowState="maximized");
subplot(3,1,1)
hold on
plot(Sow_KF,X_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,X_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,2)
hold on
plot(Sow_KF,Y_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,Y_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,3)
hold on
plot(Sow_KF,Z_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,Z_LS,"DisplayName","LS",LineStyle="--")

subplot(3,1,1)
title("ECEF坐标时间序列图")
xlabel("Sow/s")
ylabel("X/m")
legend
subplot(3,1,2)
xlabel("Sow/s")
ylabel("Y/m")
legend
subplot(3,1,3)
xlabel("Sow/s")
ylabel("Z/m")
legend

exportgraphics(img,".\pics\ECEF_Coord.png",'Resolution',600)
%% 
img = figure(WindowState="maximized");
subplot(3,1,1)
hold on
plot(Sow_KF,B_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,B_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,2)
hold on
plot(Sow_KF,L_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,L_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,3)
hold on
plot(Sow_KF,H_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,H_LS,"DisplayName","LS",LineStyle="--")

subplot(3,1,1)
title("WGS84坐标系")
xlabel("Sow/s")
ylabel("B/deg")
legend
subplot(3,1,2)
xlabel("Sow/s")
ylabel("L/deg")
legend
subplot(3,1,3)
xlabel("Sow/s")
ylabel("H/m")
legend

exportgraphics(img,".\pics\BLH_Coord.png",'Resolution',600)
%% 
img=figure(WindowState="maximized");
subplot(3,1,1)
hold on
plot(Sow_KF,E_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,E_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,2)
hold on
plot(Sow_KF,N_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,N_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,3)
hold on
plot(Sow_KF,U_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,U_LS,"DisplayName","LS",LineStyle="--")

subplot(3,1,1)
title("ENU坐标系")
xlabel("Sow/s")
ylabel("E/m")
legend
subplot(3,1,2)
xlabel("Sow/s")
ylabel("N/m")
legend
subplot(3,1,3)
xlabel("Sow/s")
ylabel("U/m")
legend

exportgraphics(img,".\pics\ENU_Coord.png",'Resolution',600)
%% 
img=figure(WindowState="maximized");
subplot(3,1,1)
hold on
plot(Sow_KF,Vx_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,Vx_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,2)
hold on
plot(Sow_KF,Vy_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,Vy_LS,"DisplayName","LS",LineStyle="--")
subplot(3,1,3)
hold on
plot(Sow_KF,Vz_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,Vz_LS,"DisplayName","LS",LineStyle="--")

subplot(3,1,1)
title("ECEF速度")
xlabel("Sow/s")
ylabel("Vx/m·s^{-1}")
legend
subplot(3,1,2)
xlabel("Sow/s")
ylabel("Vy/m·s^{-1}")
legend
subplot(3,1,3)
xlabel("Sow/s")
ylabel("Vz/m·s^{-1}")
legend

exportgraphics(img,".\pics\ECEF_Vel.png",'Resolution',600)
%% 
img=figure(WindowState="maximized");
subplot(2,1,1)
hold on
plot(Sow_KF,clock_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,clock_LS,"DisplayName","LS",LineStyle="--")
subplot(2,1,2)
hold on
plot(Sow_KF,VC_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,VC_LS,"DisplayName","LS",LineStyle="--")

subplot(2,1,1)
title("接收机钟差")
xlabel("Sow/s")
ylabel("{\Delta}Clock/m")
legend
subplot(2,1,2)
title("接收机钟速")
xlabel("Sow/s")
ylabel("V_{{\Delta}Clock}/m·s^{-1}")
legend

exportgraphics(img,".\pics\Clock.png",'Resolution',600)
%% 
img=figure(WindowState="maximized");
hold on
plot(Sow_KF,BS_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,BS_LS,"DisplayName","LS",LineStyle="--")

title("BDS卫星数量")
xlabel("Sow/s")
ylabel("Num")
legend

exportgraphics(img,".\pics\BDS_Sates_Num.png",'Resolution',600)
img=figure(WindowState="maximized");
hold on
plot(Sow_KF,GS_KF,"DisplayName","KF",LineStyle="-")
plot(Sow_LS,GS_LS,"DisplayName","LS",LineStyle="--")

title("GPS卫星数量")
xlabel("Sow/s")
ylabel("Num")
legend

exportgraphics(img,".\pics\GPS_Sates_Num.png",'Resolution',600)
%% 
img = figure(WindowState="maximized");
meanE=mean(E_KF);
meanN=mean(N_KF);
scatter(E_KF,N_KF,5,U_KF,'filled','DisplayName','E-N')
colormap('parula')
c=colorbar('Ticks',linspace(min(U_KF),max(U_KF),6));
set(get(c,'Title'),'string','U/m');
title("E-N-KF")
xlabel("E/m")
ylabel("N/m")
grid on
axis equal
hold on
scatter(meanE,meanN,40,[1,0.2,0.5],'filled','^','DisplayName','Ave-Pos')
legend
exportgraphics(img,".\pics\E-N-KF.png",'Resolution',600)

img=figure(WindowState="maximized");
meanE=mean(E_LS);
meanN=mean(N_LS);
scatter(E_LS,N_LS,5,U_LS,'filled','DisplayName','E-N')
colormap('parula')
c=colorbar('Ticks',linspace(min(U_LS),max(U_LS),6));
set(get(c,'Title'),'string','U/m');
title("E-N-LS")
xlabel("E/m")
ylabel("N/m")
grid on
axis equal
hold on
scatter(meanE,meanN,40,[1,0.2,0.5],'filled','^','DisplayName','Ave-Pos')
legend
exportgraphics(img,".\pics\E-N-LS.png",'Resolution',600)
%% 
img=figure(WindowState="maximized");
scatter(E_LS,N_LS,5,'filled','DisplayName','LS',Marker='o')
hold on
scatter(E_KF,N_KF,5,'filled','DisplayName','KF',Marker='^')
title("E-N")
xlabel("E/m")
ylabel("N/m")
grid on
axis equal
hold on
legend
exportgraphics(img,".\pics\E-N.png",'Resolution',600)
%% 
close all
