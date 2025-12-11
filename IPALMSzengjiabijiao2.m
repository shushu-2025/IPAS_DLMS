clc;        % clear command window
clear all;  % clear all variables
% p=[0,5;
%     0,0;
%     1,0;
%     1,5];


%%
N = 16;             % number of agents in the network
L = 32;%64*2;             % size of wo, which  is Mx1
wopt=randn(L,1);wopt=wopt/norm(wopt);wopt(1)=-abs(wopt(1));
load w2eq1 wopt
wopt=1*wopt;
Num_iter  = 100000;  % number of iterations per experiment
Num_trial =20%150;   % number of experiments
%load RU20;
%load A5;%A20;
A_ba=1/16*ones(N,N);
load A_lap
A_ba=A_uniform;%A_lap;%
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
MSD_lms=zeros(N,Num_iter_new);MSD_lmsno=zeros(N,Num_iter_new);
MSD_dlms=zeros(N,Num_iter_new);MSD_GT=zeros(N,Num_iter_new);
MSD_SGT=zeros(N,Num_iter_new);
MSD_MEX=zeros(N,Num_iter_new);MSD_MEXno=zeros(N,Num_iter_new);
MSD_lmsMD=zeros(N,Num_iter_new);
savez=zeros(1,Num_iter_new);saveav=zeros(1,Num_iter_new);
savez2=zeros(1,Num_iter_new);
savezno2=zeros(1,Num_iter_new);saveavno2=zeros(1,Num_iter_new);
csavecwxno2=zeros(1,Num_iter_new);
csavecw=zeros(1,Num_iter_new);csavecwx=zeros(1,Num_iter_new);
csavecwx2=zeros(1,Num_iter_new);
savezdlms=zeros(1,Num_iter_new);saveavdlms=zeros(1,Num_iter_new);
saveavdlms_np=zeros(1,Num_iter_new);
csavecwdlms=zeros(1,Num_iter_new);csavecwdlms_np=zeros(1,Num_iter_new);
%A_ba=A;%A;%0.6*eye(N)+0.4*A;
mn=0.4;var=0.0002;
bet=mn*(1-mn)^2/var-1-mn;%%v=alp*bet/((alp+bet)^2*(alp+bet+1))=m*(1-m)/(b/(1-m)+1)
alp=mn*bet/(1-mn);
mumax=1;
mu=0.0005*rand(N,1);%步长
load mu mu
mu=1/2*mu;
tau=0.8;
Rx=rand(N,1);Rv=1*rand(N,1);Rm=1*rand(N,1);
for i=1:Num_iter_new
    Rxsin(i)=(sin(2*pi*0.0002*i)+1.2)*0.1;%(sin(2*pi*0.0002*i)+1.5)*0.2;
end
load infc Rx Rv Rm
saveg=zeros(L,Num_iter_new);savemug=zeros(1,Num_iter_new);
change=Num_iter_new;
for tt=1:Num_trial % iterating over experiments
    tt
    zavsave=0;zavsave2=0;zavsavedlms=0;zavsavedlms_np=0;
    zavsaveno2=0;
    w=1*(2*rand(L,N)-1);%randn(L,N);
    wdlms=w;w2=w;
    x_vec=zeros(L,N);
    wdlms_np=w;
    wexGT=w;wno=w;w2no=w;wMD=w;
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
    pjiex=[];pjiexpro=[];pjiexpro2=[];PDLMSpjiex=[];
    for i=1:Num_iter_new
        u1(i)= sqrt(Rxsin(i))*randn;
    end
    xinput=zeros(L,1);
    for i=L:Num_iter_new
        if i <change
            for k=1:N
                %r(k)=sqrt((wx-p(k,1))^2+(wy-p(k,2))^2)+sqrt(0.001)*randn;
                v(k)=1*sqrt(1*Rv(k))*randn;
                m(k)=sqrt(1*Rm(k))*randn;
                if k==1
                    %      u(:,k)=sqrt(Rxsin(i))*randn(L,1);
                    xinput=[u1(i:-1:i-L+1)]';
                    dk(k)=v(k)+(1+m(k))*xinput'*wopt;
                    x_vec(:,k)=xinput;
                else
                    u(:,k)=sqrt(Rx(k))*randn(L,1);
                    dk(k)=v(k)+(1+m(k))*u(:,k)'*wopt;
                    x_vec(:,k)=u(:,k);
                end
            end
        end

        for k=1:N
            MSD_dlms(k,i)=MSD_dlms(k,i)+(wopt-wdlms(:,k))'*(wopt-wdlms(:,k));
            MSD_GT(k,i)=MSD_GT(k,i)+(wopt-wdlms_np(:,k))'*(wopt-wdlms_np(:,k));
            MSD_lmsno(k,i)=MSD_lmsno(k,i)+(wopt-wno(:,k))'*(wopt-wno(:,k));
            MSD_MEXno(k,i)=MSD_MEXno(k,i)+(wopt-w2no(:,k))'*(wopt-w2no(:,k));
            MSD_lmsMD(k,i)=MSD_lmsMD(k,i)+(wopt-wMD(:,k))'*(wopt-wMD(:,k));
            MSD_SGT(k,i)=MSD_SGT(k,i)+(wopt-wexGT(:,k))'*(wopt-wexGT(:,k));
        end

        % for k=1:N
        %     if k==1
        %         bex(k)=2*Rxsin(i);
        %         bex2(k)=2*Rxsin(i);
        %         jufcl(k)=2*bex(k);
        %     else
        %         bex(k)=2*Rx(k);
        %         bex2(k)=2*Rx(k);
        %         jufcl(k)=2*bex(k);
        %     end
        % end

        %% no归一化
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
                kappa=1;
            end
            muk1=mu(k)/(1+kappa*thea1(k));
            phino(:,k)=wno(:,k)-muk1*(gprim(:,k));%修改地方
        end
        for k=1:N
            noiseadd(:,k)=tau/sqrt(L)*norm(wno(:,k))*ones(L,1);
        end
        zsave=phino(1,1)+noiseadd(1,1);etaav=(i-1)/i;%0.99;%;
        savez(i)=savez(i)+zsave; zavsave=etaav*zavsave+(1-etaav)*zsave;
        csavecwx(i)=csavecwx(i)+zavsave;

        if mod(i,L)==0&i~=1
            xhat=mu(k)*g(:,1)-noiseadd(:,1)+noiseadd(:,2);xhat=xhat/(std(xhat)+1e-4);
            pjiexpro=[pjiexpro;xhat(end:-1:1)];
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
                kappa=1;
            end
            muk=mu(k)/(1+kappa*thea2(k));
            phi2no(:,k)=w2no(:,k)-muk*(gprim(:,k));%修改地方
        end
        for k=1:N
            noiseadd(:,k)=tau/sqrt(L)*norm(w2no(:,k))*ones(L,1);
        end
        zsaveno2=phi2no(1,1)+noiseadd(1,1);etaav=(i-1)/i;%0.99;%;
        savezno2(i)=savezno2(i)+zsaveno2; zavsaveno2=etaav*zavsaveno2+(1-etaav)*zsaveno2;
        csavecwxno2(i)=csavecwxno2(i)+zavsaveno2;
        if mod(i,L)==0&i~=1
            xhat=mu(k)*g2(:,1)-noiseadd(:,1)+noiseadd(:,2);xhat=xhat/(std(xhat)+1e-4);
            pjiexpro2=[pjiexpro2;xhat(end:-1:1)];
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
        zsavedlms=phidlms(1,1)+noiseadd(1,1);
        zavsavedlms=etaav*zavsavedlms+(1-etaav)*zsavedlms;%(i-1)/i;1/i
        saveavdlms(i)=saveavdlms(i)+zavsavedlms;

        csavecwdlms(i)=csavecwdlms(i)+zsavedlms(1,1);
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
        if mod(i,L)==0&i~=1
            xhat=mu(k)*gdlms(:,1)-noiseadd(:,1);xhat=xhat/(std(xhat)+1e-4);
            PDLMSpjiex=[PDLMSpjiex;xhat(end:-1:1)];
        end

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
        zsavedlms_np=phidlms_np(1,1);
        zavsavedlms_np=etaav*zavsavedlms_np+(1-etaav)*zsavedlms_np;%(i-1)/i;1/i
        saveavdlms_np(i)=saveavdlms_np(i)+zavsavedlms_np;

        csavecwdlms_np(i)=csavecwdlms_np(i)+zsavedlms_np;
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

        if mod(i,L)==0&i~=1
            xhat=mu(k)*gdlms_np(:,1);xhat=xhat/(std(xhat)+1e-4);
            pjiex=[pjiex;xhat(end:-1:1)];
        end

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
            phiMD(:,k)=wMD(:,k)-1*(gprimMD(:,k));
        end
        for k=1:N
            noiseadd(:,k)=1/sqrt(L)*norm(wMD(:,k))*ones(L,1)+1*sqrt(1*jufcl(k))*randn(L,1);
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

        %% %%%%%%%%%%%%%%%ex PP-DLMS
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
MSD_av_db_ddlms = 10*log10(MSD_av2dlms);

%%
tt2no = MSD_lmsno/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2no = sum(tt2no)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_db_dno = 10*log10(MSD_av2no);

tt2MEXno =MSD_MEXno/Num_trial; %
MSD_av2MEXno = sum(tt2MEXno)/N;   %
MSD_av_db_dMEXno = 10*log10(MSD_av2MEXno);

tt2GT =MSD_GT/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2GT = sum(tt2GT)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_db_dGT = 10*log10(MSD_av2GT);

tt2MD = MSD_lmsMD/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2MD = sum(tt2MD)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_db_dMD = 10*log10(MSD_av2MD);

tt2SGT =MSD_SGT/Num_trial; % each row contains the MSD evolution of the corresponding agent
MSD_av2SGT = sum(tt2SGT)/N;   % learning curve for adaptive combination matrix averaged over experiments and number of nodes
MSD_av_db_dSGT = 10*log10(MSD_av2SGT);
figure%提出算法1
subplot(4, 1, 3);
hold on
%plot(iter,savez,iter,saveav,'b',iter,csavecw,iter,w_s(5)*ones(1,Num_iter),'r')
%p1=plot(iter,savez,'k-');p1.MarkerSize = 10;
%p1.MarkerIndices = 1:400000:length(iter);
p2=plot(iter,savez(L:Num_iter_new-1)/Num_trial,'g-.*');p2.MarkerSize = 10;
ylim([-0.3 0.6]);
p2.MarkerIndices = 1:floor(Num_iter/6):length(iter);
p11=plot(iter,csavecwx(L:Num_iter_new-1)/Num_trial,'r-s');p11.MarkerSize = 10;
ylim([-0.3 0.6]);
p11.MarkerIndices = 1:floor(Num_iter/8):length(iter);%滑动平均
p3=plot(iter,wopt(1)*ones(1,Num_iter),'r--*');p3.MarkerSize = 10;
ylim([-0.3 0.6]);
p3.MarkerIndices = 1:floor(Num_iter/8):length(iter);
%p3=plot(iter,-wopt(1)*ones(1,Num_iter),'r--*');p3.MarkerSize = 10;
%p3.MarkerIndices = 1:300000:length(iter);
legend('${\phi}_k(n)$','sliding average of ${\phi}_{k,1}(n)$','${w}^*_{1}$');
xlabel('Iteration number ($n$)');
ylabel('Amplitude');

% plot(iter,savezdlms,iter,saveavdlms,'b',iter,csavecwdlms,iter,w_s(5)*ones(1,Num_iter),'r')
%figure %PDLMS
subplot(4, 1, 2);
hold on
%p1=plot(iter,savezdlms,'k-');p1.MarkerSize = 10;
%p1.MarkerIndices = 1:400000:length(iter);
p21=plot(iter,csavecwdlms(L:Num_iter_new-1)/Num_trial,'k-.s');p21.MarkerSize = 10;
ylim([-0.3 0.6]);
p21.MarkerIndices = 1:floor(Num_iter/6):length(iter);
p2=plot(iter,saveavdlms(L:Num_iter_new-1)/Num_trial,'g-.s');p2.MarkerSize = 12;
ylim([-0.3 0.6]);
p2.MarkerIndices = 1:floor(Num_iter/8):length(iter);
p3=plot(iter,wopt(1)*ones(1,Num_iter),'r--*');p3.MarkerSize = 10;
ylim([-0.3 0.6]);
p3.MarkerIndices = 1:floor(Num_iter/8):length(iter);
%p3=plot(iter,-wopt(1)*ones(1,Num_iter),'r--*');p3.MarkerSize = 10;
%p3.MarkerIndices = 1:300000:length(iter);
legend('${\phi}_{k,1}(n)$','sliding average of ${\phi}_{k,1}(n)$','${w}^*_{1}$');
xlabel('Iteration number ($n$)');
ylabel('Amplitude');
%figure%非保护的
subplot(4, 1, 1);
hold on
%p1=plot(iter,savezdlms,'k-');p1.MarkerSize = 10;
%p1.MarkerIndices = 1:400000:length(iter);
p21=plot(iter,csavecwdlms_np(L:Num_iter_new-1)/Num_trial,'k-.s');p21.MarkerSize = 10;
ylim([-0.5 0.6]);
p21.MarkerIndices = 1:floor(Num_iter/6):length(iter);
p2=plot(iter,saveavdlms_np(L:Num_iter_new-1)/Num_trial,'g-.s');p2.MarkerSize = 12;
ylim([-0.5 0.6]);
p2.MarkerIndices = 1:floor(Num_iter/8):length(iter);
p3=plot(iter,wopt(1)*ones(1,Num_iter),'r--*');p3.MarkerSize = 10;
ylim([-0.5 0.6]);
p3.MarkerIndices = 1:floor(Num_iter/8):length(iter);
legend('${\psi}_{k,1}(n)$','sliding average of ${\psi}_{k,1}(n)$','${w}^*_{1}$');
xlabel('Iteration number ($n$)');
ylabel('Amplitude');
%%
subplot(4, 1, 4);
hold on
%p1=plot(iter,savezdlms,'k-');p1.MarkerSize = 10;
%p1.MarkerIndices = 1:400000:length(iter);
p21=plot(iter,savezno2(L:Num_iter_new-1)/Num_trial,'k-.s');p21.MarkerSize = 10;
ylim([-0.5 0.6]);
p21.MarkerIndices = 1:floor(Num_iter/6):length(iter);
p2=plot(iter,csavecwxno2(L:Num_iter_new-1)/Num_trial,'g-.s');p2.MarkerSize = 12;
ylim([-0.5 0.6]);
p2.MarkerIndices = 1:floor(Num_iter/8):length(iter);
p3=plot(iter,wopt(1)*ones(1,Num_iter),'r--*');p3.MarkerSize = 10;
ylim([-0.5 0.6]);
p3.MarkerIndices = 1:floor(Num_iter/8):length(iter);
legend('${\phi}_{k,1}(n)$','sliding average of ${\phi}_{k,1}(n)$','${w}^*_{1}$');
xlabel('Iteration number ($n$)');
ylabel('Amplitude');

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
legend('DLMS without protection mechanism','PDLMS','MDSG','PP-DLMS','proposed IPAS-DLMS1','proposed IPAS-DLMS2');

figure
subplot(3, 2, 1:2);
plot(iter,u1(1:Num_iter_new-L))
ylim([-4 4]);xlabel('Sample number');
ylabel('Amplitude');
legend('true')
subplot(3, 2, 3);
plot(iter,pjiex(1:Num_iter_new-L))
ylim([-4 4]);xlabel('Sample number');
ylabel('Amplitude');
legend('DLMS')
subplot(3, 2, 4);
plot(iter,PDLMSpjiex(1:Num_iter_new-L))
ylim([-4 4]);xlabel('Sample number');
ylabel('Amplitude');
legend('PDLMS')
subplot(3, 2, 5);
plot(iter,pjiexpro(1:Num_iter_new-L))
ylim([-4 4]);xlabel('Sample number');
ylabel('Amplitude');
legend('proposed BCP-DLMS1')
subplot(3, 2, 6);
plot(iter,pjiexpro2(1:Num_iter_new-L))
ylim([-4 4]);xlabel('Sample number');
ylabel('Amplitude');
legend('proposed BCP-DLMS2')


DLMS=calc_mi(pjiex(1:Num_iter_new-L), u1(1:Num_iter_new-L)', 15)
PDLMS=calc_mi(PDLMSpjiex(1:Num_iter_new-L), u1(1:Num_iter_new-L)', 15)
DLMS1=calc_mi(pjiexpro(1:Num_iter_new-L), u1(1:Num_iter_new-L)', 15)
DLMS2=calc_mi(pjiexpro2(1:Num_iter_new-L), u1(1:Num_iter_new-L)', 15)


