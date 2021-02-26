%addpath /scratch/s4688360/matlab/matdcd-1.0

clear all
close all

tic;
%dirname = '/lustre/scratch2/s4688360/lmpruns/50_solvated/1/pppm/56'
% dirname = '/lustre/scratch2/s4688360/lmpruns/50_solvated/'
%dirname = '/scratch/s4688360/lmpruns/50_solvated/tc_automatic/20_10k/60more/'
%dirname = '/scratch/s4688360/lmpruns/50_solvated/tc_automatic/crosscheck/'
%dirname = '/scratch/s4688360/lmpruns/10_solvated/4bigbox10/'

%dirname = '/scratch/ws/0/s4688360-Azo/10_solvated/4bigbox10/'
dirname = './'

%dirname = '/lustre/scratch2/s4688360/lmpruns/50_solvated/switch/test/again/back'
%dirname = '/lustre/scratch2/s4688360/lmpruns/50_solvated/1/switch/'
%dirname = '/lustre/scratch2/s4688360/lmpruns/50_solvated/cl1_1_10/switch05/'
%nbins=350;
freq=500*4;
%freq=100;
%freq=25;
dt=2.0;

r_cutoff = 5.8
angle_cutoff = 15
cos_cutoff = cos(angle_cutoff*pi/180) %0.966 %0.966 %5 %66

%subdir={'/rerun_noh2o/'};
%subdir={'2/rerun12/'};
%subdir={'rerun'};
%subdir={'cluster_analysis'};
%subdir={'cluster_analysis_2nd'};
%subdir={'cluster_analysis_4th_amorphous'};

%subdir={'new_cluster_analysis_bigg'};
%subdir={'1/5.8_0.966', '2/5.8_0.966', '4/5.8_0.966', '5/5.8_0.966',  };
subdir={'example' };



%subdir={'1/rerun12/','2/rerun12/','3/rerun12/','4/rerun12/','5/rerun12/',};
%subdir={'/01ps/rr','/02ps/rr','/03ps/rr', '/05ps/rr','/1ps/rr','/2ps/rr','/3ps/rr','/5ps/rr','/10ps/rr'};
%subdir={'/01ps/back/rrback','/02ps/back/rrback','/03ps/back/rrback', '/05ps/back/rrback','/1ps/back/rrback','/2ps/back/rrback','/3ps/back/rrback','/5ps/back/rrback','/10ps/back/rrback'};
cd(dirname);

%% Verallgemeinerung der Typen fuer schnelle Anpassung

length(subdir);

%Ncom = 100;

for b=1:length(subdir)

        curfile=subdir{b}; %% current file location
	subpath=strcat(dirname,subdir{b});

        %% centers of mass (dcd file)
	
        boxfn=strcat(dirname,subdir{b},'/boxsize.txt');
        boxsz=dlmread(boxfn,'',2,1);
        
        %ellfn=strcat(dirname,subdir{b},'/ellipsoid.txt');
        ellfn=strcat(dirname,subdir{b},'/ellipsoid_tensor.txt');
        ell=dlmread(ellfn,'',[3 7 -1 12 ]);  
        ell(~any(ell,2),:)=[];
        vec=ell(:,1:3);        
        rg=ell(:,4:6);
        
        
        xyzfn=strcat(dirname,subdir{b},'/com1.out');
        xyz=dlmread(xyzfn,'',[3 1 -1  3]);  
        xyz(~any(xyz(:,2:3),2), :)=[];
        

        Ncom=size(xyz,1)/size(boxsz,1);
        Nframes=size(boxsz,1)
        %Nframes=12501;
        %Nframes=20000

        cdot=zeros(size(xyz,1)/Ncom);
        cosines=zeros(size(xyz,1)/Ncom);
        ind = [] %zeros(size(xyz,1)/Ncom);
        chist=zeros(Nframes, 100);
        ctotal=zeros(Nframes,1);
        Nas=zeros(Nframes,1);
        Nat=zeros(Nframes,1);
        Nmax=zeros(Nframes,1); % biggest cluster size at time t
        Navg=zeros(Nframes,1); % average cluster size, alternative to Nas
        
        dotprod_t = zeros(Nframes,1);
        angles_t = zeros(Nframes,1);
        dist_t = zeros(Nframes,1);
        
        all_dist = [];
        all_angles = [];
        all_dotprod = [];
        
        dotprod_avg_t = zeros(Nframes,1);
        angles_avg_t = zeros(Nframes,1);
        dist_avg_t = zeros(Nframes,1);
        
        dotprod_max_t = zeros(Nframes,1);
        dotprod_min_t = zeros(Nframes,1);
        angles_max_t = zeros(Nframes,1);
        angles_min_t = zeros(Nframes,1);        
        dist_max_t = zeros(Nframes,1);
        dist_min_t = zeros(Nframes,1);        

        
        dotprod_std_t = zeros(Nframes,1);
        angles_std_t = zeros(Nframes,1);
        dist_std_t = zeros(Nframes,1);
        
	%cind=[];

        for k=1:Nframes %50ns %1:Nframes %1:Nframes % 17000:109:17200 %1:Nframes %2001:500:20001%1:1%tsteps  %size(abc,1)/Ncom  %% welcher snapshot (zeit)

            if (mod(k, 100)==0)
                k
            end

            lx=boxsz(k,1);
            ly=boxsz(k,2);
            lz=boxsz(k,3);
            
            lxinv = 1./lx;
            lyinv = 1./ly;
            lzinv = 1./lz;
            
            indices = [];
            
            dist_collect = [];
            dotprod_collect = [];
            angles_collect = [];
            
            cind=zeros(400,2);

            clist=zeros(100,100);

            p=1;
            % double loop over all pairs
            count=0;
            
            for i=1:Ncom-1 %[23 27 64 31]%1:Ncom-1

                oivec=vec(i+(k-1)*Ncom,1:3);
                rgi=rg(i+(k-1)*Ncom,1);

                for j=i+1:Ncom %[23 27 64 31]%i+1:Ncom%Ncom
                    %[i j]

                    ojvec=vec(j+(k-1)*Ncom,1:3);       % orientation vector  
                    rgj=rg(j+(k-1)*Ncom,1:3);       

                    % distance of particles
                    % respect minimum image convention!!!!!!
                    % first fold back to box

                    %i+(k-1)*Ncom
                    %j+(k-1)*Ncom
                    
                    xi=xyz(i+(k-1)*Ncom,1);
                    xj=xyz(j+(k-1)*Ncom,1);
                    dx = abs(xj -xi);  
                    dx = dx - round(dx * lxinv) * lx;
                    
                    yi=xyz(i+(k-1)*Ncom,2);
                    yj=xyz(j+(k-1)*Ncom,2);
                    dy=abs(yj-yi);
                    dy = dy - round(dy * lyinv) * ly;

                    zi=xyz(i+(k-1)*Ncom,3);
                    zj=xyz(j+(k-1)*Ncom,3);
                    dz=zj-zi;
                    dz = dz - round(dz * lzinv) * lz;

                    r2 = dx*dx+dy*dy+dz*dz;
                    r=sqrt(r2);
                    rvec=[dx dy dz];
                    rup = rvec*oivec';
                    rupvec = rup * oivec; % /norm(oivec)
                    rsidevec = rvec - rupvec;
                    rside = norm(rsidevec);



                    dotprod=oivec * ojvec';
                    cdot(k)=dotprod;
                    cosines(k)=acos(dotprod)*180/pi;

                    %if (r < 6.0) % this works well
                    if (r < r_cutoff && dotprod > cos_cutoff)
                        
                        %if (r < 8.8)
                        count = count+1;
                        
                        [i j];
                        indices=[indices, i, j];
                        cind(p,1:2)=[i j];
                        p=p+1;
                        
                        dist_t(k) = dist_t(k) + r;
                        dotprod_t(k) = dotprod_t(k) + dotprod;
                        angles_t(k) = angles_t(k) + real(acos(dotprod)) * 180/pi;
                        
                        dist_collect = [dist_collect; r];
                        dotprod_collect = [dotprod_collect; dotprod];
                        angles_collect = [angles_collect; real(acos(dotprod)) * 180/pi];
                        
                        all_dist = [all_dist; r];
                        all_dotprod = [all_dotprod; dotprod];
                        all_angles = [all_angles; real(acos(dotprod)) * 180/pi];
                        
                        
                        
                    end


                end
                
            end
        


            %%NOW FIND CLUSTERS FROM THE PAIRS
        indices=unique(indices);
        ind=[ind,size(indices,2)];
        
        %return

        
        %return
            %% add the inverse pairs to the matrix, and sort by 1st column again (keeping paris together)
                %cind=sortrows(vertcat(cind(1:count, 1:2),[cind(1:count,2)  cind(1:count,1)]),1);
            % find the rows that contain number
        p2=1;
        ilist=1:Ncom;
        whichlist=zeros(Ncom,1);
        for i=ilist
                i;

                %return
                % check if there are more than 2 appearances
                % introduce a array of length Ncom, that contains 0 if the mol.id is not yet in a shortlist
                % or contains the number of the shortlist at the Nth position (N=the index of paired atom) of array

                if (size((find(cind==i)),1)>0 )
                    %if (size((find(cind==i)),1)==1 )
                        %return
                    %end
                    % find indices of rows (=pair-indices), where i appears
                        idx=union(find(cind(:,2)==i), find(cind(:,1)==i));
                    candidates=unique(cind(idx,:))';
                    % determine if the molindex i is already sorted into a shortlist or not
                    if (sum(whichlist(candidates))==0)
                        getindex=p2;
                        molid=unique(cind(idx,:))';
                        clist(getindex,1:max(size(molid)))=molid;
                        whichlist(molid)=p2;
                        p2=p2+1;
                    else	% CHECK IF THERE ARE DIFFERENT SHORTLISTS ALREADY MADE FOR THE ELMENTS
                        %get shortlist-id of the already considered element				
                        getindices=whichlist(candidates(find(whichlist(candidates))));
                        getindex=min(getindices);
                        molid=union(unique(cind(idx,:))', clist(getindex,:));
                        molid(molid==0)=[];
                        clist(getindex,1:max(size(molid)))=molid;
                        otherindices=getindices(getindices~=getindex);
                        clist(otherindices, :)=0;
                        whichlist(molid)=getindex;
                    end
                end
                final_list = clist(any(clist,2),:);
		%%%% now find out how many 2-cluters, 3-clusters ... we have
		cdist=zeros(1,100);

		for l=1:size(final_list,1)
			sizeidx=sum(final_list(l,:)>0);
			cdist(sizeidx)=cdist(sizeidx)+1;
		end
		chist(k,:)=cdist(:);
		ctotal(k)=sum(cdist(:));
        %numbersv=1:100;
        Nat(k)=chist(k,:)*(1:100)'; % total number of molecules that are in a cluster of any size at time t
        Nas(k)=Nat(k)/ctotal(k);    % avg cluster size (2, 3,4  molecular)
        
        if (sum(cdist==0) < 100)
            Nmax(k) = max(find(cdist));
            Navg(k) = nonzeros(chist(k, :) * (1:100)')/(sum(chist(k,:))); %mean(nonzeros(cdist .* (1:100)));
        end
        %dist_t(k) = dist_t(k)/count; % average of any pair that is somehow part of a cluster  
        %dotprod_t(k) = dotprod_t(k)/count;
        %angles_t(k) = angles_t(k)/count;
        
        end
            
        dist_t(k) = dist_t(k)/count; % average of any pair that is somehow part of a cluster  
        dotprod_t(k) = dotprod_t(k)/count;
        angles_t(k) = angles_t(k)/count;
        
        dist_avg_t(k) = mean(dist_collect);
        dist_std_t(k) = std(dist_collect);
        dotprod_avg_t(k) = mean(dotprod_collect);
        dotprod_std_t(k) = std(dotprod_collect);
        angles_avg_t(k) = mean(angles_collect);
        angles_std_t(k) = std(angles_collect);
        
        if count  > 0
            dist_min_t(k) = min(dist_collect);
            dist_max_t(k) = max(dist_collect);
            dotprod_min_t(k) = min(dotprod_collect);
            dotprod_max_t(k) = max(dotprod_collect);        
            angles_min_t(k) = min(angles_collect);
            angles_max_t(k) = max(angles_collect);        
        end
	    
        end

cd(subpath)


%% histograms


% h31=figure('Visible','on');
% %h = histogram(all_dist,' BinLimits',[0,6],'BinWidth',0.2)
% h = histogram(all_dist)
% %savefig(h,'all_dist.png', 'png')
% saveas(h31, 'all_dist_hist.png', 'png')

edge_min = 0.0;
edge_max = r_cutoff;
edge_width = 0.2;
dist_edges = [edge_min:edge_width:edge_max];

[bin_counts, edges] = histcounts(all_dist, dist_edges);
dist_total_avg = mean(all_dist)

%bin_counts
%edges
bin_centers = edges(1:end-1)+(edges(2)-edges(1))/2;

h32 = figure('Visible','off');
bar(bin_centers, bin_counts, 1.0, 'c')
hold on
plot([dist_total_avg, dist_total_avg], [0 , max(bin_counts)*1.0], 'r', 'LineWidth', 2)
saveas(h32,'dist_hist.png', 'png')
dlmwrite('dist_hist.txt', [bin_centers; bin_counts]', 'delimiter', '\t',  'precision', 8);


edge_min = 0.0;
edge_max = angle_cutoff;
edge_width = 0.5; %0.5;
dist_edges = [edge_min:edge_width:edge_max];

[bin_counts, edges] = histcounts(all_angles, dist_edges);
angles_total_avg = mean(all_angles)
bin_centers = edges(1:end-1)+(edges(2)-edges(1))/2;

h33 = figure('Visible','off');
bar(bin_centers, bin_counts, 1.0, 'c')
hold on
plot([angles_total_avg, angles_total_avg], [0 , max(bin_counts)*1.0], 'r', 'LineWidth', 2)
saveas(h33,'angles_hist.png', 'png')
dlmwrite('angles_hist.txt', [bin_centers; bin_counts]', 'delimiter', '\t',  'precision', 8);

%return  


%% PLOT AND SAVE AS IMAGES

steps2ns=freq/1000000 * dt;

h1=figure('Visible','off');
plot((0:Nframes-1).*steps2ns, ctotal(:), 'k', (0:Nframes-1).*steps2ns, chist(:,2), 'r', (0:Nframes-1).*steps2ns, chist(:,3), 'b', (0:Nframes-1).*steps2ns, chist(:,4), 'g',  (0:Nframes-1).*steps2ns, chist(:,5), 'm', (0:Nframes-1).*steps2ns, chist(:,6), 'c', (0:Nframes-1).*steps2ns, chist(:,7), 'y', 'Linewidth',1);
title('Cluster size distribution over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Number of clusters of a specific size','Fontsize',16);
legend('# of clusters total','Clusters of 2', 'of 3', 'of 4', 'of 5', 'of 6', 'of 7', 'Location','northwest');
legend boxoff;
saveas(h1,'ctotal.png','png');
saveas(h1,'ctotal.eps','epsc');
dlmwrite('ctotal.txt', [0:Nframes-1; ctotal'; chist(:,2)'; chist(:,3)'; chist(:,4)'; chist(:,5)'; chist(:,6)'; chist(:,7)'; chist(:,8)'; chist(:,9)'; chist(:,10)';]', 'delimiter', '\t',  'precision', 8);

h2=figure('Visible','off'); %figure(2)
plot((0:Nframes-1).*steps2ns, Nat(:),'Linewidth',1);
title('Number of atoms in all clusters over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Number of atoms in all clusters','Fontsize',16);
saveas(h2,'Nat.png','png');
saveas(h2,'Nat.eps','epsc');
dlmwrite('Nat.txt', [0:Nframes-1; Nat']', 'delimiter', '\t',  'precision', 8);

h12=figure('Visible','off'); %figure(2)
plot((0:Nframes-1).*steps2ns, Nmax(:),'Linewidth',1);
title('Largest Cluster Size over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Largest cluster','Fontsize',16);
saveas(h12,'Nmax.png','png');
saveas(h12,'Nmax.eps','epsc');
dlmwrite('Nmax.txt', [0:Nframes-1; Nmax']', 'delimiter', '\t',  'precision', 8);


h3=figure('Visible','off'); %figure(3)
plot((0:Nframes-1).*steps2ns, Nas(:), 'Linewidth',1);
title('Average Cluster Size over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Average Cluster Size','Fontsize',16);
axis([0 (Nframes-1)*steps2ns 0 max(Nas)+1])
saveas(h3,'Nas.png','png');
saveas(h3,'Nas.eps','epsc');
dlmwrite('Nas.txt', [0:Nframes-1; Nas']', 'delimiter', '\t',  'precision', 8);

h13=figure('Visible','off'); %figure(3)
plot((0:Nframes-1).*steps2ns, Navg(:), 'Linewidth',1);
title('Average Cluster Size over time Altern.','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Average Cluster Size Alt.','Fontsize',16);
axis([0 (Nframes-1)*steps2ns 0 max(Navg)+1])
saveas(h13,'Navg.png','png');
saveas(h13,'Navg.eps','epsc');
dlmwrite('Navg.txt', [0:Nframes-1; Nas']', 'delimiter', '\t',  'precision', 8);


h7=figure('Visible','off'); %figure(3)
%plot((0:Nframes-1).*steps2ns, dist_avg_t(:), 'b', 'Linewidth',1);
plot((0:Nframes-1).*steps2ns, dist_min_t(:), 'r', 'Linewidth',1);
hold on;
plot((0:Nframes-1).*steps2ns, dist_max_t(:), 'g', 'Linewidth',1);
hold on;
plot((0:Nframes-1).*steps2ns, dist_avg_t(:), 'b', 'Linewidth',1);
title('Average Distance between Pairs','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Average Pair Distance','Fontsize',16);
axis([0 (Nframes-1)*steps2ns 0 max(dist_avg_t)*1.15])
saveas(h7,'dist_t.png','png');
saveas(h7,'dist_t.eps','epsc');
dlmwrite('dist_t.txt', [0:Nframes-1; dist_avg_t'; dist_std_t']', 'delimiter', '\t',  'precision', 8);

h8=figure('Visible','off'); %figure(3)
%plot((0:Nframes-1).*steps2ns, dotprod_avg_t(:), 'b' ,'Linewidth',1);
plot((0:Nframes-1).*steps2ns, dotprod_min_t(:), 'r','Linewidth',1);
hold on;
plot((0:Nframes-1).*steps2ns, dotprod_max_t(:), 'g','Linewidth',1);
hold on;
plot((0:Nframes-1).*steps2ns, dotprod_avg_t(:), 'b' ,'Linewidth',1);
title('Average Cosine between Pairs','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Average Pair Cosine','Fontsize',16);
axis([0 (Nframes-1)*steps2ns min(dotprod_min_t(:)) 1.1])
saveas(h8,'dotprod_t.png','png');
saveas(h8,'dotprod_t.eps','epsc');
dlmwrite('dotprod_t.txt', [0:Nframes-1; dotprod_avg_t'; dotprod_std_t']', 'delimiter', '\t',  'precision', 8);

h9=figure('Visible','off'); %figure(3)
%plot((0:Nframes-1).*steps2ns, angles_avg_t(:), 'b' , 'Linewidth',1);
plot((0:Nframes-1).*steps2ns, angles_min_t(:), 'r','Linewidth',1);
hold on;
plot((0:Nframes-1).*steps2ns, angles_max_t(:), 'g','Linewidth',1);
hold on;
plot((0:Nframes-1).*steps2ns, angles_avg_t(:), 'b' , 'Linewidth',1);
title('Average Angle between Pairs','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Average Pair Angle','Fontsize',16);
axis([0 (Nframes-1)*steps2ns 0 max(angles_avg_t)*1.2])
saveas(h3,'Nas.png','png');
saveas(h3,'Nas.eps','epsc');
dlmwrite('angles_t.txt', [0:Nframes-1; angles_avg_t'; angles_std_t']', 'delimiter', '\t',  'precision', 8);




%return


%return

%avgwindow=freq;
%avgwindow=freq/(10*2);
avgwindow=freq/(10*4);
sat=size(Nat,1);
Nat_avg=zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);
Nat_crd=zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);
crd=zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);
Nas_avg=zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);
ctotal_avg=zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);
chist_avg=zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),100);
Nmax_avg = zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);

dist_window =  zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);
dotprod_window =  zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);
angles_window =  zeros(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+~(mod(sat,avgwindow)==1),1);


for ii=1:fix(sat/avgwindow)
	chist_avg(ii,:)=mean(chist(avgwindow*(ii-1)+1:avgwindow*ii,:));
	crd(ii)=mean((avgwindow*(ii-1)+1:avgwindow*ii));
	Nat_avg(ii)=mean(Nat(avgwindow*(ii-1)+1:avgwindow*ii));
	Nat_crd(ii)=mean((avgwindow*(ii-1)+1:avgwindow*ii));
	Nas_avg(ii)=mean(Nas(avgwindow*(ii-1)+1:avgwindow*ii));
	ctotal_avg(ii)=mean(ctotal(avgwindow*(ii-1)+1:avgwindow*ii));
    Nmax_avg(ii)=mean(Nmax(avgwindow*(ii-1)+1:avgwindow*ii));
    
    dist_window(ii)=mean(dist_avg_t(avgwindow*(ii-1)+1:avgwindow*ii));
    dotprod_window(ii)=mean(dotprod_avg_t(avgwindow*(ii-1)+1:avgwindow*ii));
    angles_window(ii)=mean(angles_avg_t(avgwindow*(ii-1)+1:avgwindow*ii));
	
end
if (~(~mod(sat,avgwindow)))==1
	
	if (mod(sat,avgwindow))==1
		chist_avg(fix(sat/avgwindow)+1,:)=chist(end,:);
	else
		chist_avg(fix(sat/avgwindow)+(~(~mod(sat,avgwindow))),:)=mean(chist(avgwindow*fix(sat/avgwindow)+1:end,:));
		chist_avg(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1,:)=chist(end,:);	
	end
	Nat_avg(fix(sat/avgwindow)+1)=mean(Nat(avgwindow*fix(sat/avgwindow)+1:end));
	Nat_crd(fix(sat/avgwindow)+1)=mean((avgwindow*fix(sat/avgwindow)+1:size(Nat,1)));
	Nas_avg(fix(sat/avgwindow)+1)=mean(Nas(avgwindow*fix(sat/avgwindow)+1:end));
	crd(fix(sat/avgwindow)+1)=mean((avgwindow*fix(sat/avgwindow)+1:size(Nat,1)));
	ctotal_avg(fix(sat/avgwindow)+1)=mean(ctotal(avgwindow*fix(sat/avgwindow)+1:end));
    Nmax_avg(fix(sat/avgwindow)+1)=mean(Nmax(avgwindow*fix(sat/avgwindow)+1:end));
    
    dist_window(fix(sat/avgwindow)+1)=mean(dist_avg_t(avgwindow*fix(sat/avgwindow)+1:end));
    dotprod_window(fix(sat/avgwindow)+1)=mean(dotprod_avg_t(avgwindow*fix(sat/avgwindow)+1:end));
    angles_window(fix(sat/avgwindow)+1)=mean(angles_avg_t(avgwindow*fix(sat/avgwindow)+1:end));
    
end
if ~(mod(sat,avgwindow)==1)
	Nat_avg(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=Nat(end);
	Nas_avg(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=Nas(end);
	crd(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=sat;
	Nat_crd(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=sat;
	ctotal_avg(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=ctotal(end);
    Nmax_avg(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=Nmax(end);
    
    dist_window(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=dist_avg_t(end);
    dotprod_window(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=dotprod_avg_t(end);
    angles_window(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=angles_avg_t(end);
end


Nat_avg(:) = round(Nat_avg(:));
Nas_avg(:) = round(Nas_avg(:));
ctotal_avg(:) = round(ctotal_avg(:));
chist_avg = round(chist_avg);

Nmax_avg = round(Nmax_avg);
%Nat_avg(fix(sat/avgwindow)+(~(~mod(sat,avgwindow)))+1)=Nat(end);




%chist_avg=[chist(1,:) chist_avg(:,:)'];
chist_avg=[chist(1,:); chist_avg(:,:)];
Nat_avg=[Nat(1) Nat_avg(:)'];
Nat_crd=[0 Nat_crd(:)'];
crd=[0 crd(:)'];
Nas_avg=[Nas(1) Nas_avg(:)'];
Nmax_avg=[Nmax(1) Nmax_avg(:)'];

dist_window=[dist_avg_t(1) dist_window(:)'];
dotprod_window=[dotprod_avg_t(1) dotprod_window(:)'];
angles_window=[angles_avg_t(1) angles_window(:)'];

ctotal_avg=[ctotal(1) ctotal_avg(:)'];
%dlmwrite('chist_avg.txt', [crd; chist_avg]', 'delimiter', '\t',  'precision', 8);





h4=figure('Visible','off');
plot(crd.*steps2ns, ctotal_avg(:), 'k', crd.*steps2ns, chist_avg(:,2), 'r', crd.*steps2ns, chist_avg(:,3), 'b', crd.*steps2ns, chist_avg(:,4), 'g',  crd.*steps2ns, chist_avg(:,5), 'm', crd.*steps2ns, chist_avg(:,6), 'c', crd.*steps2ns, chist_avg(:,7), 'y', 'Linewidth',1);
title('Cluster size distribution over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Number of clusters of a specific size','Fontsize',16);
legend('# of clusters total','Clusters of 2', 'of 3', 'of 4', 'of 5', 'of 6', 'of 7', 'Location','northwest');
legend boxoff;
xlim([0 crd(end)*steps2ns])
saveas(h4,'ctotal_avg.png','png');
saveas(h4,'ctotal_avg.eps','epsc');

h5=figure('Visible','off'); %figure(2)
plot(crd.*steps2ns, Nat_avg(:),'Linewidth',1);
title('Number of atoms in all clusters over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
xlim([0 crd(end)*steps2ns])
ylabel('Number of atoms in all clusters','Fontsize',16);
saveas(h5,'Nat_avg.png','png');
saveas(h5,'Nat_avg.eps','epsc');

%h6=figure('Visible','off'); %figure(3)
h6=figure('Visible','off'); %figure(3)
plot(crd.*steps2ns, Nas_avg(:), 'Linewidth',1);
title('Average Cluster Size over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Average Cluster Size','Fontsize',16);
axis([0 crd(end)*steps2ns 0 5.0])
saveas(h6,'Nas_avg.png','png');
saveas(h6,'Nas_avg.eps','epsc');

%h6=figure('Visible','off'); %figure(3)
h16=figure('Visible','off'); %figure(3)
plot(crd.*steps2ns, Nmax_avg(:), 'Linewidth',1);
title('Max Cluster Size over time','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('Max Cluster Size','Fontsize',16);
axis([0 crd(end)*steps2ns 0 max(Nmax_avg)+1])
saveas(h16,'Nmax_avg.png','png');
saveas(h16,'Nmax_avg.eps','epsc');


h18=figure('Visible','off'); %figure(3)
plot(crd.*steps2ns, dist_window(:), 'Linewidth',1);
title('avg pair distances (rolling avg)','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('avg pair dist (rolling)','Fontsize',16);
axis([0 crd(end)*steps2ns 0 max(dist_window)*1.2])
saveas(h18,'dist_avg_avg.png','png');
saveas(h18,'dist_avg_avg.eps','epsc');

h19=figure('Visible','off'); %figure(3)
plot(crd.*steps2ns, dotprod_window(:), 'Linewidth',1);
title('avg pair cosine (rolling avg)','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('avg pair cosine (rolling)','Fontsize',16);
axis([0 crd(end)*steps2ns 0 max(dotprod_window)*1.2])
saveas(h19,'dotprod_avg_avg.png','png');
saveas(h19,'dotprod_avg_avg.eps','epsc');

h20=figure('Visible','off'); %figure(3)
plot(crd.*steps2ns, angles_window(:), 'Linewidth',1);
title('avg pair angles (rolling avg)','Fontsize',16);
xlabel('t [ns]','Fontsize',16);
ylabel('avg pair angles (rolling)','Fontsize',16);
axis([0 crd(end)*steps2ns 0 max(angles_window)*1.2])
saveas(h20,'angles_avg_avg.png','png');
saveas(h20,'angles_avg_avg.eps','epsc');




%return

%fname = fullfile(subpath, 'chist_avg.txt');
%fname
%fid=fopen(fname, 'wt');
%fprintf(fid, '%s\n','#Time evolution of number of clusters made of 2, 3, ... molecules'); %header
%fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', '#step','t [ns]', 'cl-of-2','3','4', '5', '6', '7', '8', '9', '10');  % header
%fclose(fid);
dlmwrite('chist_avg.txt', [crd;crd.*steps2ns ;chist_avg(:,2:10)']', 'delimiter', '\t',  'precision', 4, '-append');

%fname = fullfile(subpath, 'Nat_avg.txt');
%fid=fopen(fname, 'wt');
%fprintf(fid, '%s\n','#Number of molecules that are in a cluster at time t'); %header
%fprintf(fid, '%s\t%s\t%s\t%s\n','#step','t [ns]', 'N', 'N/Ntotal'); %header
%fclose(fid);
dlmwrite('Nat_avg.txt', [Nat_crd;Nat_crd.*steps2ns ; Nat_avg; Nat_avg./Ncom]', 'delimiter', '\t',  'precision', 8, '-append');

%fname = fullfile(subpath, 'Nas_avg.txt');
%fid=fopen(fname, 'wt');
%fprintf(fid, '%s\n','#Average cluster size (# of molecules per cluster) at time t'); %header
%fprintf(fid, '%s\t%s\t%s\n','#step','t [ns]', 'N_avg'); %header
%fclose(fid);
dlmwrite('Nas_avg.txt', [crd; crd.*steps2ns ; Nas_avg]', 'delimiter', '\t',  'precision', 8, '-append');

%fname = fullfile(subpath, 'ctotal_avg.txt');
%fid=fopen(fname, 'wt');
%fprintf(fid, '%s\n','#Number of clusters (of any size) at time t'); %header
%fprintf(fid, '%s\t%s\t%s\n','#step','t [ns]', 'N_cl'); %header
%fclose(fid);
dlmwrite('ctotal_avg.txt', [crd;crd.*steps2ns ; ctotal_avg]', 'delimiter', '\t',  'precision', 8, '-append');


dlmwrite('Nmax_avg.txt', [crd;crd.*steps2ns ; Nmax_avg]', 'delimiter', '\t',  'precision', 8, '-append');


dlmwrite('dist_avg_avg.txt', [crd;crd.*steps2ns ; dist_window]', 'delimiter', '\t',  'precision', 8, '-append');
dlmwrite('dotprod_avg_avg.txt', [crd;crd.*steps2ns ; dotprod_window]', 'delimiter', '\t',  'precision', 8, '-append');
dlmwrite('angles_avg_avg.txt', [crd;crd.*steps2ns ; angles_window]', 'delimiter', '\t',  'precision', 8, '-append');

%figure(5)
%plot(Nat_crd, Nat_avg,'g','Linewidth',4);

end
%figure(2)
   % plot(1:20001, ctotal(k), 'r', 1:20001, chist(:,2), 'b', 1:20001, chist(k,3), 'g', 1:20001, chist(k,4), 'm');
toc
