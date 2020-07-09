%% Martensite Variants
%
%%
% In this section we discuss the austenite (fcc) to ferrite (bcc) phase transformation on the example of 
% an EBSD data set collected on a plessitic microstructure of the Emsland iron meteorite.
% Plessite is the greek description for filling iron and occurs as remaining volumes between the already
% transformed kamacite (bcc in meteorites) rims. Plessite regionons are commonly surrounded by a very thin taenite (fcc) ribbons. 
% The filling iron contains as major phases again bcc and fcc, where the orientation of fcc practically always 
% indicates the orientation of the formerly huge fcc grain in the planetary body which can easily reach the dimension of meters.

plotx2east

% import the ebsd data
mtexdata emsland

% extract crystal symmetries
cs_bcc = ebsd('Fe').CS;
cs_aus = ebsd('Aus').CS;

% recover grains
ebsd = ebsd('indexed');

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);
ebsd(grains(grains.grainSize<=2)) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);

grains = smooth(grains,4);

%%
% The following lines plot bcc according to the crystallographic description
% of the selected reference direction (IPF coloring), whereas austeniteis displayed as phase in blue.

plot(ebsd('Fe'),ebsd('Fe').orientations)
hold on
plot(grains.boundary,'lineWidth',2,'lineColor','gray')
plot(grains('Aus'),'FaceColor','blue','edgeColor','b','lineWidth',1,'DisplayName','Austenite')
hold off

%%
% As expected, we recognize very small remaining fcc grains. 
% This high-temperatur phase is stabilized by the increasing nickel content 
% during transformation. The low-temperature bcc phase can solve in maximum only 
% 6\% nickel so that fcc has to assimilate the excess nickel. Size and amount of fcc is therefore
% and indication of the overall nickel content. 
% Considering only the parent fcc phase and plotting the
% orientations into an axis angle plot

plot(ebsd('Aus').orientations,'axisAngle')

%%
% we recognize the uniform orientation of all fcc grains. Deviations are assumed to be the result 
% of deformations during high-speed collisions in asteroitic belt.
% We can get this parent grain orientation by taking the 
% <orientation.mean.html |mean|> and compute the fit by the command
% <orientation.std.html |std|>

parenOri = mean(ebsd('Aus').orientations)

fit = std(ebsd('Aus').orientations) ./ degree

%%
% Next we display the bcc orientations (blue dots) in pole figures, and additionally we plot on
% top of them the parent fcc orientation (red dots).

childOri = grains('Fe').meanOrientation;

h_bcc = Miller({1,0,0},{1,1,0},{1,1,1},cs_bcc);
h_fcc = Miller({1,0,0},{1,1,0},{1,1,1},cs_aus);

plotPDF(childOri,h_bcc,'MarkerSize',5);

nextAxis(1)
hold on
plot(parenOri * h_fcc(1).symmetrise ,'MarkerFaceColor','r')
xlabel('$(100)$','Color','red','Interpreter','latex')

nextAxis(2)
plot(parenOri * h_fcc(3).symmetrise ,'MarkerFaceColor','r')
xlabel('$(111)$','Color','red','Interpreter','latex')

nextAxis(3)
plot(parenOri * h_fcc(2).symmetrise ,'MarkerFaceColor','r')
xlabel('$(110)$','Color','red','Interpreter','latex')
hold off

drawNow(gcm)

%%
% The partial coincidence of bcc and fcc poles suggests 
% an existing crystallographic orientation relationship between both phases. 
% For this technically relevant phase transformation several models have been developed 
% in order to realize what might happen on an atomic scale. One of the first 
% assumption was the Kurdjumov-Sachs (KS) orientation relationship model. It hypothesizes 
% a transition of one of the {111}-fcc into one of the {110}-bcc 
% planes since these are the closed-packed planes and the atomic arrangement is very similar. 
% Moreover, within these planes one of the <110> directions 
% of fcc is assumed to remain parallel to one of the <111> directions
% of the formed bcc. 
% Fortunately, for cubic crystals identically indexed (hkl) and [uvw]
% generate the same directions. Therefore, the pole figures for lattice planes look identical
% to those of equally indexed lattice directions, so that the same pole figure can be used 
% for the evaluation of lattice directions as well as planes.
%
% Although we could alternatively apply the MTEX command 
% |orientation.KurdjumovSachs(cs_aus,cs_bcc)|, let us define the upper mentioned orientation
% relationship manually:

KS = orientation.map(Miller(1,1,1,cs_aus),Miller(0,1,1,cs_bcc),...
      Miller(-1,0,1,cs_aus),Miller(-1,-1,1,cs_bcc));

% Please note the condition of perpendicular vectors for both phases: hu+kv+lw=0, i.e., the indices need to
% be correctly defined.

plotPDF(variants(KS,parenOri),'add2all','MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',2)

%%
% As common for models, they never match exactly the experiment. Whereas often 
% statistical reasons prevent a better match, the orientation relationship very likely depends on the 
% chemical composition at the interface. So it is no miracle that KS does not perfectly fit to the experimental
% data.
% For a quantification of the deviation between the model and the discovered conditions in plessite, 
% as simplest indicator we can first compute the mean angular deviation between all parent-to-child
% misorientaitons and the KS model.

% Each parent-to-child misorientations can be calculated by
mori = inv(childOri) * parenOri;

% whereas the mean angular deviation (output in degree) results from
mean(angle(mori, KS)) ./ degree

%fit = sqrt(mean(min(angle_outer(childOri,variants(KS,parenOri)),[],2).^2))./degree


%% Estimating the parent to child orientation relationship
%
% Since KS shows a certain deviation, we may have asked ourselfs whether there is an orientation relationship
% that better matches the measured data. 
% A canocial candidate would be the <orientation.mean.html |mean|> of all
% misorientations.

% The mean of all measured parent-to-child misorientations
p2cMean = mean(mori,'robust')

plotPDF(childOri,h_bcc,'MarkerSize',5);
hold on
plotPDF(variants(p2cMean,parenOri),'add2all','MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',2)
hold off

% mean angular deviation in degree
mean(angle(mori, p2cMean)) ./ degree

%%
% Here we have made use of our comfortable situation to know the parent orientation. 
% However, if this reference orientation is unknown we may still estimate
% the parent-to-child orientation relationship soleley from the child-to-child 
% misorientations by an algorithm derived from Tuomo Nyyssönen and available in MTEX
% by the function <calcParent2Child.html |calcParent2Child|>. This
% iterative algorithm needs as starting point an orientation relationship
% close to the actual one. Instead of KS, we can use the Nishiyama-Wassermann (NW)
% orientation relationship which also uses the closed-packed planes but not the close-packed
% directions. Instead, a symmetrical alignement is preferred which assumes one of three <110> 
% in fcc parallel to the <001> in bcc. This reduces the amount of possibilities (variants) per 
% lattice plane to three compared to six for KS. Multiplied by the number of symmetric {111} planes
% for KS result 24 variants, whereas for NW only 12 are possible. Despite the different crystallographic 
% indexing, the misorientation angle between both models is below 6 degrees.


% MTEX already contains a definition of NW which can be simply called:
NW = orientation.NishiyamaWassermann(cs_aus,cs_bcc);

% With the following command  all child-to-child misorientations will be extracted: 
grainPairs = neighbors(grains('Fe'));
ori = grains(grainPairs).meanOrientation;

% The estimation of the parent-to-child orientation relationship results from
p2cIter = calcParent2Child(ori,NW)

% and the mean angular deviation is
mean(angle(mori,p2cIter)) ./degree

%%
% The parent-to-child orientation relationship solely computed
% from the child-to-child misorientations actually fits the experimental
% orientation relationship quite well. 
%
%% Classification of variants
%
% Once we have determined parent orientations and an 
% orientation relationship we may proceed further by classifying the 
% different variants. This is computed by the command
% <calcChildVariant.html |calcChildVariant|>.

% determine for each child orientation a variantId
[variantId, packetId] = calcChildVariant(parenOri,childOri,p2cIter);

% colorize the orientations according to the variantID
color = ind2color(variantId);
plotPDF(childOri,color,h_bcc,'MarkerSize',5);

%%
% While it is very hard to distinguish the different variants in the pole
% figure plots it becomes more clear in an axis angle plot

plot(childOri,color,'axisAngle')

%%
% Another classification of variants is a grouping of them
% according to one of the assumed transformation planes {111}. 

color = ind2color(packetId);
plotPDF(childOri,color,h_bcc,'MarkerSize',5,'points',1000);

%%
% Another classification of variants is a grouping of them
% according to one of the assumed transformation planes {111}. 
% For colors - red, blue, orange and green - are sufficient to generate  
% a pole distribution wich represents each a group of variants (packet).

nextAxis(1)
hold on
opt = {'MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',3};
plot(parenOri * h_fcc(1).symmetrise ,opt{:})
xlabel('$(100)$','Color','red','Interpreter','latex')

nextAxis(2)
plot(parenOri * h_fcc(3).symmetrise ,opt{:})
xlabel('$(111)$','Color','red','Interpreter','latex')

nextAxis(3)
plot(parenOri * h_fcc(2).symmetrise ,opt{:})
xlabel('$(110)$','Color','red','Interpreter','latex')
hold off

drawNow(gcm)

%%
% We may also use the packet colors to vizualize the spatial distribution 
% of them in the EBSD map.

plot(grains('Fe'),color)

