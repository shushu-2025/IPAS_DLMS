clc;        % clear command window
clear all;  % clear all variables
% p=[0,5;
%     0,0;
%     1,0;
%     1,5];


%%
N = 16;             % number of agents in the network
L = 32;%64*2;             % size of wo, which is Mx1
wopt=randn(L,1);wopt=wopt/norm(wopt);
load w2eq1 wopt
Num_iter  = 5000000;  % number of iterations per experiment
Num_trial =1;%150;   % number of experiments
%load RU20;
%load A5;%A20;
A_ba=1/16*ones(N,N);
load A_lap
A_ba=A_lap;%A_uniform;%
[V, D] = eig(A_ba);       % 对转置矩阵求特征向量（左特征向量=转置矩阵右特征向量的转置）
eig_vals = diag(D);        % 提取特征值
[~, idx] = min(abs(eig_vals - 1));  % 定位最接近1的特征值索引

% 3. 提取特征向量、去负、归一化（分量和为1）
vx = V(:, idx);             % 原始特征向量（列向量）
vx = abs(vx);                % 取绝对值确保正分量（Perron向量性质）
vx = vx / sum(vx);            % 归一化（概率向量）
p = vx';                    % 转为行向量（左特征向量形式）

%A_ba=[1/2,1/4, 0, 1/4;
%    1/4,1/2, 1/4, 0;
%    0,  1/4, 1/2, 1/4;
%    1/4, 0,  1/4, 1/2];
%A_ba=0.99*A_ba+0.01*eye(N);
%load sigma_out20;
%load wo_m5;
%load w_s5;
% wo_m=randn(L,N);
% load wo_m5yigou
% sumRu=0;sumRuW=0;
% c=1;
% for k=1:N
%     sumRu=sumRu+RU(1,k)*eye(L);
%     sumRuW=sumRuW+RU(1,k)*eye(L)*c*wo_m(:,k);
% end
% w_s=inv(sumRu)*sumRuW;w_s(1)=-0.46;
Num_iter_new=Num_iter+L;
% mu=0.0007;%0.000001
jufcl=(0.5)*ones(1,N);%rand(1,N);
MSD_lmsno=zeros(N,Num_iter_new);
MSD_dlms=zeros(N,Num_iter_new);MSD_GT=zeros(N,Num_iter_new);
MSD_MEXno=zeros(N,Num_iter_new);
MSD_SGT=zeros(N,Num_iter_new);MSD_lmsMD=zeros(N,Num_iter_new);


%A_ba=A;%A;%0.6*eye(N)+0.4*A;
mn=0.4;var=0.0002;
bet=mn*(1-mn)^2/var-1-mn;%%v=alp*bet/((alp+bet)^2*(alp+bet+1))=m*(1-m)/(b/(1-m)+1)
alp=mn*bet/(1-mn);
mumax=1;
mu=0.001*rand(N,1);%步长
load mu mu
mu=0.0005./(N*p);
tau=0.8;
Rx=rand(N,1);Rv=1*rand(N,1);Rm=1*rand(N,1);
for i=1:Num_iter_new+L
    Rxsin(i)=(sin(2*pi*0.0002*i)+1.2)*0.1;
end
load infc Rx Rv Rm
saveg=zeros(L,Num_iter_new+L);savemug=zeros(1,Num_iter_new+L);
change=Num_iter_new+L;
for tt=1:Num_trial % iterating over experiments
    tt
    zavsave=0;zavsave2=0;zavsavedlms=0;zavsavedlms_np=0;
    w=1*(2*rand(L,N)-1);%randn(L,N);
    wdlms=w;w2=w;
    x_vec=zeros(L,N);
    wdlms_np=w;wMD=w;
    wexGT=w;wno=w;w2no=w;
    gg=zeros(L,N);
    deltFp=zeros(L,N);
    
    wmex=w;phimexbef=w;
    bex=0.5*ones(N,1);bex2=0.5*ones(N,1);
    thea2=0*ones(N,1);thea1=0*ones(N,1);
    thea1no=0*ones(N,1);thea2no=0*ones(N,1);
    dk=zeros(1,N);avtd1=0*ones(L,N);
    avetd2=0*ones(N,1);avetd1=0*ones(N,1);
    avetd1no=0*ones(N,1); avetd2no=0*ones(N,1);
    avtd2=0*ones(L,N);saveinpu=zeros(Num_iter_new,1);
    pjiex=[];
    for i=1:Num_iter_new
        u1(i)= sqrt(Rxsin(i))*randn;
    end
     xinput=zeros(L,1);
    for i=L:Num_iter_new
        if i <change
            for k=1:N
                %r(k)=sqrt((wx-p(k,1))^2+(wy-p(k,2))^2)+sqrt(0.001)*randn;
                v(k)=1*sqrt(Rv(k))*randn;
                m(k)=sqrt(Rm(k))*randn;
                if k==1
                    %                     u(:,k)=sqrt(Rxsin(i))*randn(L,1);;xinput(1:end-1)
                    xinput=[u1(i:-1:i-L+1)]';xinput(end)=xinput(end-1);%1;%
                    dk(k)=v(k)+(1+m(k))*xinput'*wopt;
                    x_vec(:,k)=xinput;
                else
                    u(:,k)=sqrt(Rx(k))*randn(L,1);u(end,k)=u(end-1,k);%1;%
                    dk(k)=v(k)+(1+m(k))*u(:,k)'*wopt;
                    x_vec(:,k)=u(:,k);
                end
                
                %r(k)=u(:,k)'*(wopt-p(k,:)')+sqrt(0.001)*randn;
                
                %%
                %                 uo(:,k)=(wopt'-p(k,:))/norm(wopt'-p(k,:),2);%arccos(uo(:,k)'*[1;0]/norm(uo(:,k)))+pi/2
                %                 u(:,k)=uo(:,k)+1*sqrt(0.1)*randn(1)*(R*uo(:,k))+1*sqrt(0.01)*randn(1)*uo(:,k);
                %                 dk(k)=u(:,k)'*wopt+1*sqrt(0.5)*randn;
                
            end
            %saveinpu(i)=u(1,1);
            %         else
            %             for k=1:N
            %                 x_vec(:,k)=sqrt(RU(1,k)).*randn(L,1);
            %                 y =  -c*w_s'*x_vec(:,k);
            %                 %var_y=sum(y.^2)/Num_iter_new;
            %                 sigma_vo(k)=1/(10^(sigma_out(k)/10));
            %                 dk(k)=y+randn(1,1)*sqrt(sigma_vo(k));
            %             end
        end
        for k=1:N
            if k==1
                mr(:,:,k)=diag(Rxsin(i:-1:i-L+1)); 
                mr(L,L-1,k)=Rxsin(i-L+2);mr(L-1,L,k)=Rxsin(i-L+2);mr(L,L,k)=Rxsin(i-L+2);
            else
                mr(:,:,k)=diag([Rx(k)*ones(L,1)]); mr(L,L-1,k)=Rx(k);mr(L-1,L,k)=Rx(k);
            end
        end
        pxt=kron(p,eye(L));
        for k=1:N
            MSD_dlms(k,i)=MSD_dlms(k,i)+(wopt-sum(pxt*reshape(wdlms,L*N,1),2)/1)'*mr(:,:,k)*(wopt-sum(pxt*reshape(wdlms,L*N,1),2)/1);
            MSD_GT(k,i)=MSD_GT(k,i)+(wopt-sum(pxt*reshape(wdlms_np,L*N,1),2)/1)'*mr(:,:,k)*(wopt-sum(pxt*reshape(wdlms_np,L*N,1),2)/1);
            MSD_lmsno(k,i)=MSD_lmsno(k,i)+(wopt-sum(pxt*reshape(wno,L*N,1),2)/1)'*mr(:,:,k)*(wopt-sum(pxt*reshape(wno,L*N,1),2)/1);
            MSD_MEXno(k,i)=MSD_MEXno(k,i)+(wopt-sum(pxt*reshape(w2no,L*N,1),2)/1)'*mr(:,:,k)*(wopt-sum(pxt*reshape(w2no,L*N,1),2)/1);
            MSD_lmsMD(k,i)=MSD_lmsMD(k,i)+(wopt-sum(pxt*reshape(wMD,L*N,1),2)/1)'*mr(:,:,k)*(wopt-sum(pxt*reshape(wMD,L*N,1),2)/1);
            MSD_SGT(k,i)=MSD_SGT(k,i)+(wopt-sum(pxt*reshape(wexGT,L*N,1),2)/1)'*mr(:,:,k)*(wopt-sum(pxt*reshape(wexGT,L*N,1),2)/1);
        end
        
        %% no
        for k=1:N
            x_hzvecno(:,k)=x_vec(:,k)+sqrt(bex(k))*randn(L,1);
            g(:,k)=-(x_hzvecno(:,k)*(dk(k)-x_hzvecno(:,k)'*wno(:,k))+(bex(k)*eye(L))*wno(:,k));%-(RU(1,k)*eye(L))*(c*w_s-w(:,k));%
            noi=0*randn(L,1);
            gprim(:,k)=g(:,k)+noi;
            thea1no(k)=0.99*thea1no(k)+0.01*norm(wno(:,k))^2;
            %avtd1(:,k)= 0.99*avtd1(:,k)+0.01*(x_vec(:,k)*(dk(k)-x_vec(:,k)'*w(:,k)));
            avetd1no(k)= 0.999*avetd1no(k)+0.001*(dk(k)-x_vec(:,k)'*wno(:,k));
            if norm(avetd1no(k))^2>0*0.004
                kappa=0;
            else
                kappa=0;
            end
            muk1=mu(k)/(1+kappa*thea1(k));
            phino(:,k)=wno(:,k)-muk1*(gprim(:,k));%修改地方
        end
        for k=1:N
            noiseadd(:,k)=tau/sqrt(L)*norm(wno(:,k))*ones(L,1);
        end
        for k=1:N
            tempno=0;
            for ll=1:N
                if A_ba(ll,k)~=0
                    tempno=tempno+A_ba(ll,k)*(phino(:,ll)+noiseadd(:,ll));
                end
            end
            w_tmpno(:,k)=tempno-noiseadd(:,k);
        end
        wno=w_tmpno;
        %% 2 no
        for k=1:N
            c_hzvec2no(:,k)=sqrt(bex2(k))*randn(L,1);
            g2(:,k)=-(x_vec(:,k)*(dk(k)-x_vec(:,k)'*w2no(:,k))+(dk(k)-x_vec(:,k)'*w2no(:,k))*c_hzvec2no(:,k));%-(RU(1,k)*eye(L))*(c*w_s-w(:,k));%
            gprim(:,k)=g2(:,k);
            thea2no(k)=0.99*thea2no(k)+0.01*norm(w2no(:,k))^2;
            %avtd2(:,k)= 0.995*avtd2(:,k)+0.005*(x_vec(:,k)*(dk(k)-x_vec(:,k)'*w2(:,k)));
            avetd2no(k)= 0.999*avetd2no(k)+0.001*(dk(k)-x_vec(:,k)'*w2no(:,k));
            if norm(avetd2no(k))^2>0*0.0002
                kappa=0;
            else
                kappa=0;
            end
            muk=mu(k)/(1+kappa*thea2(k));
            phi2no(:,k)=w2no(:,k)-muk*(gprim(:,k));%修改地方
        end
        for k=1:N
            noiseadd(:,k)=tau/sqrt(L)*norm(w2no(:,k))*ones(L,1);
        end
        for k=1:N
            temp2no=0;
            for ll=1:N
                if A_ba(ll,k)~=0
                    temp2no=temp2no+A_ba(ll,k)*(phi2no(:,ll)+noiseadd(:,ll));
                end
            end
            w_tmp2no(:,k)=temp2no-noiseadd(:,k);
        end
        w2no=w_tmp2no;
        
        %% PDLMS
        for k=1:N
            if i <change
                gdlms(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wdlms(:,k));%-(RU(1,k)*eye(L))*(c*w_s-wdlms(:,k));%
            else
                gdlms(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wdlms(:,k));%-(RU(1,k)*eye(L))*(-c*w_s-wdlms(:,k));%
            end
            %mu=unifrnd(0,1,L,1)*mumax;
            %mu=betarnd(alp, bet, L, 1)*0.1;%
            phidlms(:,k)=wdlms(:,k)-mu(k)*gdlms(:,k);
        end
        for k=1:N
            noiseadd(:,k)=1*sqrt(jufcl(k))*randn(L,1);
        end
        for k=1:N
            tempdlms=0;
            for ll=1:N
                if A_ba(ll,k)~=0
                    tempdlms=tempdlms+A_ba(ll,k)*(phidlms(:,ll)+noiseadd(:,ll));
                end
            end
            w_tmpdlms(:,k)=tempdlms;
        end
        wdlms=w_tmpdlms;
        %% DLMS 非保护
        for k=1:N
            if i <change
                gdlms_np(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wdlms_np(:,k));%-(RU(1,k)*eye(L))*(c*w_s-wdlms(:,k));%
            else
                gdlms_np(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wdlms_np(:,k));%-(RU(1,k)*eye(L))*(-c*w_s-wdlms(:,k));%
            end
            %mu=unifrnd(0,1,L,1)*mumax;
            %mu=betarnd(alp, bet, L, 1)*0.1;%
            phidlms_np(:,k)=wdlms_np(:,k)-mu(k)*gdlms_np(:,k);
        end
        for k=1:N
            tempdlms_np=0;
            for ll=1:N
                if A_ba(ll,k)~=0
                    tempdlms_np=tempdlms_np+A_ba(ll,k)*(phidlms_np(:,ll));
                end
            end
            w_tmpdlms_np(:,k)=tempdlms_np;
        end
        wdlms_np=w_tmpdlms_np;
        
         %% MD
        for k=1:N
            %x_vec(:,k)=[xk(i,k);x_vec(1:end-1,k)];
            if i <change
                gMD(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wMD(:,k));%-(RU(1,k)*eye(L))*(c*w_s-wMD(:,k));%
            else
                gMD(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wMD(:,k));%-(RU(1,k)*eye(L))*(-c*w_s-wMD(:,k));%
            end
            %%%mu=betarnd(alp, bet, L, 1)*0.1;%betarnd(alpha, beta, sampleSize, 1);
            %mean=alpha/(alpha+beta),Var=alpha*beta/((alpha+beta)^2*(alpha+beta+1))
            %alp=mn*bet/(1-mn);be=mn*(1-mn)^2/var-1-mn;%%v=alp*bet/((alp+bet)^2*(alp+bet+1))=m*(1-m)/(b/(1-m)+1)
            gprimMD(:,k)=(sqrt(0.01)*randn(L,1)+1/0.01*mu(k)).*gMD(:,k);%6e-05
            phiMD(:,k)=wMD(:,k)-1/1*(gprimMD(:,k));
        end
        for k=1:N
            noiseadd(:,k)=tau/sqrt(L)*norm(wMD(:,k))*ones(L,1)+1*sqrt(1*jufcl(k))*randn(L,1);
        end
        %zsave=phiMD(1,1)+noiseadd(1,1);
        %savez(i)=zsave;       % zavsave=0.9995*zavsave+0.0005*zsave;saveav(i)=zavsave;
        %csavecw(i)=phiMD(1,1);csavecwx(i)=wMD(1,1);
        for k=1:N
            temp=0;
            for ll=1:N
                if A_ba(ll,k)~=0
                    temp=temp+A_ba(ll,k)*(phiMD(:,ll)+noiseadd(:,ll));
                end
            end
            w_tmpMD(:,k)=temp-noiseadd(:,k);
        end
        wMD=(1-0.01)*wMD+0.01*w_tmpMD;
        
        %% %%%%%%%%%%%%%%%ex P
        for k=1:N
            if i <change
                deltF(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wexGT(:,k));%-(RU(1,k)*eye(L))*(c*wo_m(:,k)-wexGT(:,k));%
            else
                deltF(:,k)=-x_vec(:,k)*(dk(k)-x_vec(:,k)'*wexGT(:,k));%-(RU(1,k)*eye(L))*(-c*wo_m(:,k)-wexGT(:,k));%
            end
        end
        for k=1:N
            wexjGT(:,k)=wexGT(:,k)-mu(k)*deltF(:,k);%5e-05
        end
  
        for k=1:N               
            loc= randperm(L, 12);
            chos=zeros(L,1);
            chos(loc)=1;choshouy(:,k)=chos;
            Chos=diag(chos);
            shach(:,k) =Chos*(wexjGT(:,k)+1*sqrt(1*jufcl(k))*randn(L,1));
        end
        cphi_exGT=zeros(L,N);
        for k=1:N
            for ell=1:N        
                if A_ba(ell,k)~=0%&&(ell~=k)
                    cphi_exGT(:,k)=cphi_exGT(:,k)+A_ba(ell,k)*(shach(:,ell)...
                        +(eye(L)-diag(choshouy(:,ell)))*wexjGT(:,k));%
%                 elseif A_ba(ell,k)~=0&&ell==k
%                     cphi_exGT(:,k)=cphi_exGT(:,k)+A_ba(ell,k)*;
                end
            end
        end
        wexGT=cphi_exGT;
        
    end
end
iter = 1:Num_iter;

%%
tt2dlms = MSD_dlms/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2dlms = sum(tt2dlms)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av2dlmsz=0;
for ii=L:Num_iter_new
    MSD_av2dlmsz=MSD_av2dlmsz*(ii-1)/ii+MSD_av2dlms(ii)/ii;
    MSD_av2dlmszz(ii)=MSD_av2dlmsz;
end
MSD_av_db_ddlms = 10*log10(MSD_av2dlmszz);

%%
tt2no = MSD_lmsno/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2no = sum(tt2no)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av2noz=0;
for ii=L:Num_iter_new
    MSD_av2noz=MSD_av2noz*(ii-1)/ii+MSD_av2no(ii)/ii;
    MSD_av2nozz(ii)=MSD_av2noz;
end
MSD_av_db_dno = 10*log10(MSD_av2nozz);

tt2MEXno =MSD_MEXno/Num_trial; %
MSD_av2MEXno = sum(tt2MEXno)/N;   %
MSD_av2MEXnoz=0;
for ii=L:Num_iter_new
    MSD_av2MEXnoz=MSD_av2MEXnoz*(ii-1)/ii+MSD_av2MEXno(ii)/ii;
    MSD_av2MEXnozz(ii)=MSD_av2MEXnoz;
end
MSD_av_db_dMEXno = 10*log10(MSD_av2MEXnozz);

tt2GT =MSD_GT/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2GT = sum(tt2GT)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av2GTz=0;
for ii=L:Num_iter_new
    MSD_av2GTz=MSD_av2GTz*(ii-1)/ii+MSD_av2GT(ii)/ii;
    MSD_av2GTzz(ii)=MSD_av2GTz;
end
MSD_av_db_dGT = 10*log10(MSD_av2GTzz);%非保护DLMS

tt2MD = MSD_lmsMD/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2MD = sum(tt2MD)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av2GTMD=0;
for ii=L:Num_iter_new
    MSD_av2GTMD=MSD_av2GTMD*(ii-1)/ii+MSD_av2MD(ii)/ii;
    MSD_av2GTMDm(ii)=MSD_av2GTMD;
end
MSD_av_db_dMD = 10*log10(MSD_av2GTMDm);

tt2SGT =MSD_SGT/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2SGT = sum(tt2SGT)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av2SGTm=0;
for ii=L:Num_iter_new
    MSD_av2SGTm=MSD_av2SGTm*(ii-1)/ii+MSD_av2SGT(ii)/ii;
    MSD_av2GSTMDm(ii)=MSD_av2SGTm;
end
MSD_av_db_dSGT = 10*log10(MSD_av2GSTMDm);


figure;hold on
p112=plot(iter,MSD_av_db_dGT(L:Num_iter_new-1),'g-.s');p112.MarkerSize = 10;
p112.MarkerIndices = 1:Num_iter/10:length(iter);%非保护DLMS

p111=plot(iter,MSD_av_db_ddlms(L:Num_iter_new-1),'k-.s');p111.MarkerSize = 10;
p111.MarkerIndices = 1:Num_iter/10:length(iter);%保护DLMS

p111x=plot(iter,MSD_av_db_dMD(L:Num_iter_new-1),'g-.*');p111x.MarkerSize = 10;
p111x.MarkerIndices = 1:Num_iter/10:length(iter);

p111x1=plot(iter,MSD_av_db_dSGT(L:Num_iter_new-1),'g-.*');p111x1.MarkerSize = 10;
p111x1.MarkerIndices = 1:Num_iter/10:length(iter);

p115=plot(iter,MSD_av_db_dno(L:Num_iter_new-1),'r-.x');p115.MarkerSize = 10;
p115.MarkerIndices = 1:Num_iter/10:length(iter);
p116=plot(iter,MSD_av_db_dMEXno(L:Num_iter_new-1),'b-.x');p116.MarkerSize = 10;
p116.MarkerIndices = 1:floor(Num_iter/8):length(iter);



%plot(iter,MSD_av_db_ddlms,iter,MSD_av_db_dGT,iter,MSD_av_db_d,iter,MSD_av_db_dMEX);
xlabel('Iteration number ($n$)');
ylabel('MSD (dB)');
legend('DLMS','PDLMS','MDSG','PP-DLMS','Proposed IPAS-DLMS1','Proposed IPAS-DLMS2');


%save resu20251209




