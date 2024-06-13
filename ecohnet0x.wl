(* ::Package:: *)

BeginPackage["ecohnet0x`"]

(*LIMITS*)
LinearRegression::usage=""
GenerateXY::usage=""
GenerateXYn::usage=""
GenerateDdbDtg::usage=""
LIMITS::usage=""
LIMall::usage=""
LIMallpara::usage=""
(*CCM*)
Onestepah::usage=""
EmbQ::usage=""
Ccm::usage=""
Ccmrsmp::usage=""
CCMkendall::usage=""
CCMks::usage=""
CCMebisu::usage=""
CCMseason::usage=""
CCMallparasrg::usage=""
CCMallparanonsrg::usage=""
Seasonf::usage=""
DumFT::usage=""
EBZ::usage=""
DumSeason::usage=""
(*PCM*)
Falseneighbour::usage=""
FEmbQ::usage=""
Mcm::usage=""
Pcm::usage=""
Pcc::usage=""
SimplexProjection::usage=""
Extcon::usage=""
Mcm2Pcm::usage=""
Mcm2Pcmp::usage=""
PCMall::usage=""
PCMallpara::usage=""
(*Correkation*)
Corf::usage=""
(*EcohNet*)
SetVariables::usage=""
RCsetup::usage=""
ES::usage=""
Prop::usage=""
Reservoirep::usage=""
Reservoirepd::usage=""
BRMSE::usage=""
RMSE::usage=""
Evp::usage=""
(*misc*)
rndROCx::usage=""
std::usage=""
Norstd::usage=""
Nordiffstd::usage=""
Logstd::usage=""
Rnkstd::usage=""
LogRnkstd::usage=""
Rnkmat::usage=""
Symm::usage=""
Nev::usage=""

(*LIMITS*)
Unprotect[LinearRegression];
Unprotect[GenerateXY];
Unprotect[GenerateXYn];
Unprotect[GenerateDdbDtg];
Unprotect[LIMITS];
Unprotect[LIMall];
Unprotect[LIMallpara];
(*CCM*)
Unprotect[Onestepah];
Unprotect[EmbQ];
Unprotect[Ccm];
Unprotect[Ccmrsmp];
Unprotect[CCMkendall];
Unprotect[CCMks];
Unprotect[CCMebisu];
Unprotect[CCMseason];
Unprotect[CCMallparasrg];
Unprotect[CCMallparanonsrg];
Unprotect[Seasonf]
Unprotect[DumFT];
Unprotect[EBZ];
Unprotect[DumSeason];
(*PCM*)
Unprotect[Falseneighbour];
Unprotect[FEmbQ];
Unprotect[Mcm];
Unprotect[Pcm];
Unprotect[Pcc];
Unprotect[SimplexProjection];
Unprotect[Extcon];
Unprotect[Mcm2Pcm];
Unprotect[Mcm2Pcmp];
Unprotect[PCMall];
Unprotect[PCMallpara];
(*Correkation*)
Unprotect[Corf];
(*EcohNet*)
Unprotect[SetVariables];
Unprotect[RCsetup];
Unprotect[ES];
Unprotect[Prop];
Unprotect[Reservoirep];
Unprotect[Reservoirepd];
Unprotect[BRMSE];
Unprotect[RMSE];
Unprotect[Evp];
Unprotect[rcparameters];
Unprotect[p1];
Unprotect[p2];
Unprotect[p3];
Unprotect[p4];
Unprotect[q1];
Unprotect[q2];
Unprotect[InitializeNN];
Unprotect[comp];
Unprotect[compd];
Unprotect[RCcore];
Unprotect[RCprd];
Unprotect[RCmat];
Unprotect[RCbaselineprop];
Unprotect[RCcoren];
Unprotect[RCprdn];
Unprotect[RCmatn];
Unprotect[RCbaselinepropn];
(*misc*)
Unprotect[rndROCx];
Unprotect[std];
Unprotect[Norstd];
Unprotect[Nordiffstd];
Unprotect[Logstd];
Unprotect[Rnkstd];
Unprotect[LogRnkstd];
Unprotect[Rnkmat];
Unprotect[Symm];
Unprotect[Nev];

Begin["ecohnet0x`Private`"]


(*methods*)

(*LIMITS*)
(*functions*)
LinearRegression[iactive_,i_,data_]:=Module[{activePositions,xdb,ydb,xtg,ytg,
cs,yStar},
activePositions=Flatten[Position[iactive,1]];
{xdb,ydb,xtg,ytg}={data[[1,All,activePositions]],data[[2,All,i]],
data[[3,All,activePositions]],data[[4,All,i]]};
cs=PseudoInverse[xdb] . ydb;
yStar=(cs . #)&/@(xtg);
{N[Total[(yStar-ytg)^2]/Length[ytg]](*MSE*),
ReplacePart[iactive,MapThread[(#1->#2)&,{activePositions,cs}]](*cij^**)}]

GenerateXY[timeseries_]:=Module[{mean,xs,ys},
mean=Mean[Exp[timeseries]];
xs=Drop[Transpose[Transpose[Exp[timeseries]]-mean],-1];
ys=Drop[timeseries,1]-Drop[timeseries,-1];
{xs,ys}]

GenerateXYn[timeseries_]:=Module[{mean,xs,ys},
mean=Mean[timeseries];
xs=Drop[Transpose[Transpose[timeseries]-mean],-1];
ys=Drop[timeseries,1]-Drop[timeseries,-1];
{xs,ys}]

GenerateDdbDtg[data_]:=Module[{dbLength,xs,ys,timeStamps,
fordb,dbtk,rndStp,xdb,ydb,fortg,xtg,ytg},
{xs,ys}=data;
timeStamps=Range[Length[xs]];
dbLength=IntegerPart[Length[timeStamps]*0.5]+1;
dbtk=dbLength;
rndStp=RandomSample[timeStamps];
fordb=Take[rndStp,dbtk];
xdb=xs[[fordb]];
ydb=ys[[fordb]];
fortg=Drop[rndStp,dbtk];
xtg=xs[[fortg]];
ytg=ys[[fortg]];
{xdb,ydb,xtg,ytg}]

LIMITS[timeseries_,i_,qc_,baggingIteration_]:=Module[{baggingResult,numberofSpecies,
ddbDtg,iactive,initialResult,msePrev,mseBest,ciStar,
links,q,newIactives,results,mses,hBest},
numberofSpecies=Length[timeseries[[1]]];
baggingResult=Table[
ddbDtg=GenerateDdbDtg[GenerateXY[timeseries]];
iactive=ReplacePart[Table[0,{numberofSpecies}],i->1];
initialResult=LinearRegression[iactive,i,ddbDtg];
msePrev=initialResult[[1]];
ciStar=initialResult[[2]];
links=1;q=99;
While[links<numberofSpecies-1&&q>qc,
newIactives=ReplacePart[iactive,#->1]&/@Flatten[Position[iactive,0]];
results=(LinearRegression[#,i,ddbDtg])&/@newIactives;
mses=results[[All,1]];
mseBest=Min[mses];
q=(1-mseBest/msePrev);
hBest=Position[mses,mseBest][[1,1]];
If[q>qc,iactive=newIactives[[hBest]];msePrev=mseBest;ciStar=results[[hBest,2]]];
links++;
];
ciStar,{baggingIteration}];
{Median[baggingResult],
1-N[(Count[#,_?(#!=0&)]/Length[#])&/@Transpose[baggingResult]]}
]

LIMall[timeseries_,th_:0.]:=Module[{mm},
Transpose[Table[LIMITS[Log[timeseries],w,th,256][[{1,2}]],
{w,1,Length[timeseries[[1]]]}]]]

LIMallpara[timeseries_,th_:0.]:=Module[{mm},
Transpose[ParallelTable[LIMITS[Log[timeseries],w,th,256][[{1,2}]],
{w,1,Length[timeseries[[1]]]},DistributedContexts->Automatic]]]

(*CCM*)

(*Onestep ahead prediction*)
Onestepah[timeseries_]:=Table[Ccm[Drop[timeseries,-1],Drop[timeseries,1],w],{w,1,12}];

(*embedding dim*)
EmbQ[timeseries_]:=Module[{fn,mi},
fn=Abs[Onestepah[timeseries]];
mi=Max[fn];
Position[fn,_?(#>mi*0.95&)][[1,1]]];

(*single ccm*)
Ccm[samplex_,sampley_,embx_]:=Module[{db,yys,dbact,xt,nrm,ptc,ui,wi},
db=Transpose[Partition[#,embx,1]&/@{samplex,sampley}];
yys=Table[
dbact=Drop[db,{t}];
xt=db[[t,1]];
nrm=Norm/@Transpose[Transpose[dbact[[All,1]]]-xt];
ptc=Take[Sort[Transpose[{nrm,dbact[[All,2]]}]],embx+1];
ui=Exp[(-(ptc[[All,1]]/(ptc[[1,1]]/.(0.->1))))/._?(#<-20&)->-20];
wi=ui/(Total[ui]/.(0.->1));
{db[[t,2,-1]],wi . ptc[[All,2,-1]]},
{t,1,Length[db]}];
(Correlation[yys[[All,1]],yys[[All,2]]]+RandomReal[{-1,1}]*10^-8)
]

(*ccm core w resampling*)
Ccmrsmp[db_,sn_]:=Module[{embx=Length[db[[1,1]]],yys,dbact,xt,nrm,ptc,ui,wi},
yys=Table[
dbact=RandomSample[Drop[db,{t}],sn];
xt=db[[t,1]];
nrm=Norm/@Transpose[Transpose[dbact[[All,1]]]-xt];
ptc=Take[Sort[Transpose[{nrm,dbact[[All,2]]}]],embx+1];
ui=N[Exp[(-ptc[[All,1]]/(ptc[[1,1]]/.(0.->1)))/._?(#<-20&)->-20]];
wi=ui/(Total[ui]/.(0.->1));
{db[[t,2,-1]],wi . ptc[[All,2,-1]]},
{t,1,Length[db]}];
(Correlation[yys[[All,1]],yys[[All,2]]]+RandomReal[{-1,1}]*10^-8)
]

(*CCM with Ebisuzaki*)
CCMebisu[samplex_,sampley_,embx_,th_]:=Module[{db,ccm,qccm,nboot,minccmsmp,maxccmsmp,qmax,p1,p2,
surrogsmp,rdb,qmaxs,p3,maxccm,surrogccm},
db=Transpose[Partition[#,embx,1]&/@{samplex,sampley}];
nboot=500;
minccmsmp=Table[Ccmrsmp[db,embx+1],{nboot}];
maxccm=Ccm[samplex,sampley,embx];
surrogccm=Table[Ccm[EBZ[samplex],EBZ[sampley],embx],{nboot}];
qmaxs=Quantile[surrogccm,0.95];
p1=(1-Count[minccmsmp,_?(#<=maxccm&)]/nboot);
p2=(1-Count[surrogccm,_?(#<=maxccm&)]/nboot);
{If[Max[p1,p2]<=th,maxccm,0],N[Max[p1,p2]]}
]

(*CCM with seasonal surrogate*)
CCMseason[samplex_,sampley_,embx_,th_,cyc_]:=Module[{db,ccm,qccm,nboot,minccmsmp,maxccmsmp,qmax,p1,p2,
surrogsmp,rdb,qmaxs,p3,maxccm,surrogccm},
db=Transpose[Partition[#,embx,1]&/@{samplex,sampley}];
nboot=500;
minccmsmp=Table[Ccmrsmp[db,embx+1],{nboot}];
maxccm=Ccm[samplex,sampley,embx];
surrogccm=Table[Ccm[DumSeason[samplex,cyc],DumSeason[sampley,cyc],embx],{nboot}];
qmaxs=Quantile[surrogccm,0.95];
p1=(1-Count[minccmsmp,_?(#<=maxccm&)]/nboot);
p2=(1-Count[surrogccm,_?(#<=maxccm&)]/nboot);
{If[Max[p1,p2]<=th,maxccm,0],N[Max[p1,p2]]}
]

(*CCM with kendall*)
CCMkendall[samplex_,sampley_,embx_,th_]:=Module[{db,ccm,qccm,nboot,minccmsmp,maxccmsmp,qmax,p1,p2,
q25ccmsmp,q50ccmsmp,rdb,qmaxs,p3,maxccm},
db=Transpose[Partition[#,embx,1]&/@{samplex,sampley}];
nboot=500;
minccmsmp=Table[Ccmrsmp[db,embx+1],{nboot}];
q25ccmsmp=Table[Ccmrsmp[db,IntegerPart[Length[db]*0.25]],{nboot}];
q50ccmsmp=Table[Ccmrsmp[db,IntegerPart[Length[db]*0.5]],{nboot}];
maxccmsmp=Table[Ccmrsmp[db,IntegerPart[Length[db]*0.75]],{nboot}];
maxccm=Ccm[samplex,sampley,embx];
p1=1-(Count[minccmsmp,_?(#<=maxccm&)]/nboot);
p2=KendallTauTest[{1,2,3,4,5},
Join[Quantile[#,0.5]&/@{minccmsmp,q25ccmsmp,q50ccmsmp,maxccmsmp},{maxccm}]];
{If[Max[p1,p2]<=th,maxccm,0],N[Max[p1,p2]]}
]

(*CCM with kendall + season*)
CCMks[samplex_,sampley_,embx_,th_,cyc_]:=Module[{db,ccm,qccm,nboot,minccmsmp,maxccmsmp,qmax,p1,p2,
q25ccmsmp,q50ccmsmp,rdb,qmaxs,p3,maxccm},
db=Transpose[Partition[#,embx,1]&/@{Seasonf[samplex,cyc],Seasonf[sampley,cyc]}];
nboot=500;
minccmsmp=Table[Ccmrsmp[db,embx+1],{nboot}];
q25ccmsmp=Table[Ccmrsmp[db,IntegerPart[Length[db]*0.25]],{nboot}];
q50ccmsmp=Table[Ccmrsmp[db,IntegerPart[Length[db]*0.5]],{nboot}];
maxccmsmp=Table[Ccmrsmp[db,IntegerPart[Length[db]*0.75]],{nboot}];
maxccm=Ccm[samplex,sampley,embx];
p1=1-(Count[minccmsmp,_?(#<=maxccm&)]/nboot);
p2=KendallTauTest[{1,2,3,4,5},
Join[Quantile[#,0.5]&/@{minccmsmp,q25ccmsmp,q50ccmsmp,maxccmsmp},{maxccm}]];
{If[Max[p1,p2]<=th,maxccm,0],N[Max[p1,p2]]}
]

CCMallparasrg[timeseries_,th_,season_,cyc_:1]:=Module[{embxs,ccmm,xy},
embxs=EmbQ[#]&/@Transpose[timeseries];xy=Flatten[Table[Table[If[x==y,0,{x,y}],{x,1,Length[embxs]}],{y,1,Length[embxs]}],1];
Which[season==1,(*CCM with surrogate(Seasonality)*)
ccmm=Partition[ParallelMap[(If[Length[#]==0,{0,0},CCMseason[timeseries[[All,#[[1]]]],timeseries[[All,#[[2]]]],embxs[[#[[1]]]],th,cyc]])&,xy,DistributedContexts->Automatic],
Length[embxs]],
True,(*CCM with surrogate(Ebisuzaki)*)
ccmm=Partition[ParallelMap[(If[Length[#]==0,{0,0},CCMebisu[timeseries[[All,#[[1]]]],timeseries[[All,#[[2]]]],embxs[[#[[1]]]],th]])&,xy,DistributedContexts->Automatic],
Length[embxs]]];
{ccmm[[All,All,1]]+IdentityMatrix[Length[timeseries[[1]]]],
(ccmm[[All,All,2]])}
]

CCMallparanonsrg[timeseries_,th_,season_,cyc_:1]:=Module[{embxs,ccmm,xy},
embxs=EmbQ[#]&/@Transpose[timeseries];xy=Flatten[Table[Table[If[x==y,0,{x,y}],{x,1,Length[embxs]}],{y,1,Length[embxs]}],1];
Which[season==1,(*CCM with Kendall tau + sasonality)*)
ccmm=Partition[ParallelMap[(If[Length[#]==0,{0,0},CCMks[timeseries[[All,#[[1]]]],timeseries[[All,#[[2]]]],embxs[[#[[1]]]],th,cyc]])&,xy,DistributedContexts->Automatic],
Length[embxs]],
True,(*CCM with kendall*)
ccmm=Partition[ParallelMap[(If[Length[#]==0,{0,0},CCMkendall[timeseries[[All,#[[1]]]],timeseries[[All,#[[2]]]],embxs[[#[[1]]]],th]])&,xy,DistributedContexts->Automatic],
Length[embxs]]];
{ccmm[[All,All,1]]+IdentityMatrix[Length[timeseries[[1]]]],
(ccmm[[All,All,2]])}
]

(*Surrogate*)

DumFT[data_]:=Module[{ll,xxx,absak,aaa},(ll=Length[data];
If[OddQ[ll],xxx=Drop[data,1];ll=Length[xxx],xxx=data;ll=Length[xxx]];
absak=Abs[Fourier[xxx,FourierParameters->{-1, 1}]];
2*(ll^0.5)*Re[InverseFourier[
aaa=ReplacePart[Take[absak,ll/2]*Exp[I*RandomReal[{0,2*Pi},ll/2]],
{1->0,ll/2->(2^0.5)*absak[[ll/2]]*Cos[RandomReal[{0,2*Pi}]]}];
Join[aaa,Table[0,{ll/2}]]
]])]

EBZ[xx_]:=Module[{ip,sim,dumsim},
ip=Interpolation[DeleteCases[MapIndexed[({#2[[1]],#1})&,xx],_?(!NumericQ[#[[2]]]&)],
Method->"Spline",InterpolationOrder->1];
sim=ip[Range[Length[xx]]];
dumsim=DumFT[sim];
(If[#[[1]]=="NA","NA",#[[2]]])&/@Transpose[{xx,dumsim}]]

DumSeason[xx_,cyc_]:=Module[{mid,srep,rmid},mid=Mod[Range[Length[xx]],cyc];
srep=(#[[1,1]]->Mean[#[[All,2]]])&/@Split[Sort[Transpose[{mid,xx}]],(#[[1]]==#2[[1]]&)];
rmid=mid/.srep;
rmid+RandomSample[xx-rmid]]

Seasonf[da_,cyc_]:=Flatten[(#-Mean[#])&/@Partition[da,cyc]]

(*********************************************************************************************************)

(*PCM//Leng et al *)

(*False neighbour method*)
Falseneighbour[timeseries_]:=Module[{db,result,dx,ss},
Quiet[result=Table[db=Partition[Partition[timeseries,ww,1],2,1];
Mean[DeleteCases[(dx=#;
ss=Drop[Sort[{Norm[#[[1]]-dx[[1]]],#}&/@db],1][[1]];
(Norm[db[[1,2]]-ss[[2,2]]])/(Norm[db[[1,1]]-ss[[2,1]]]+10^-8))&/@db,Indeterminate]],
{ww,2,12}]];
result/.(Indeterminate->Mean[DeleteCases[result,Indeterminate]])
];

(*embedding dim*)
FEmbQ[timeseries_]:=Module[{fn,mi},
fn=Falseneighbour[timeseries];
mi=Min[fn];
Position[fn,_?(#<mi*1.05&)][[1,1]]];

Mcm[samplex_,sampley_,embx_]:=Module[{dbxy,y,yox},
dbxy=Transpose[Partition[#,embx,1]&/@{samplex,sampley}];(*db>explanatory, responce*)
{y,yox}=SimplexProjection[dbxy];
Correlation[y,yox]
]

Pcm[samplexyz_,embxz_]:=Module[{
samplex,sampley,samplez,embx,embz,dbxy,dbxz,y,yox,z,zox,
dbzoxy,yy,yozox,le},
{samplex,sampley,samplez}=samplexyz;
{embx,embz}=embxz;
dbxy=Transpose[Partition[#,embx,1]&/@{samplex,sampley}];(*db>explanatory, responce*)
dbxz=Transpose[Partition[#,embx,1]&/@{samplex,samplez}];
{y,yox}=SimplexProjection[dbxy];
{z,zox}=SimplexProjection[dbxz];
dbzoxy=Transpose[Partition[#,embz,1]&/@{zox,Take[sampley,-Length[zox]]}];
{yy,yozox}=SimplexProjection[dbzoxy];
le=Length[yy];
N[(Correlation[yy,Take[yox,-le]]-Correlation[yy,yozox]*Correlation[Take[yox,-le],yozox]
)/(((1-Correlation[yy,yozox]^2)*(1-Correlation[Take[yox,-le],yozox]^2))^0.5)]
]

Pcc[pxyz_]:=Module[{px,py,pz},
{px,py,pz}=pxyz/.{_?(#>1&)->0.99,_?(#<-1&)->-0.99};
N[(px-(py*pz))/(((1-((py)^2))*(1-((pz)^2)))^0.5)]/.{_?(#>1&)->0.99,_?(#<-1&)->-0.99}
]

SimplexProjection[db_]:=Module[{embx,dbp,xt,nrm,ptc,ui,wi},Transpose[Table[
embx=Length[db[[1,1]]];(*explanatory, responce*)
dbp=Drop[db,{t}];
xt=db[[t,1]];
nrm=Norm/@Transpose[Transpose[dbp[[All,1]]]-xt];
ptc=Take[Sort[Transpose[{nrm,dbp[[All,2]]}]],embx+1];
ui=Exp[(-(ptc[[All,1]]/(ptc[[1,1]]/.(0.->1))))/._?(#<-20&)->-20];
wi=ui/(Total[ui]/.(0.->1));
{db[[t,2,-1]],wi . ptc[[All,2,-1]]},
{t,1,Length[db]}]]]

Extcon[con_]:=Module[{crit,ncon,ae,be,ce,nncon},
ncon=con;
crit=Length[ncon[[1]]];
While[crit>3,
nncon=Flatten[(ae=Take[#,2];
be=Take[#,{3,-2}];
ce=Take[#,-1];
{Join[ae,be],
Join[ae[[{1}]],ce,be],
Join[ae[[{2}]],ce,be]}
)&/@ncon,1];
ncon=nncon;
crit=Length[ncon[[1]]]];
ncon]

Mcm2Pcm[mcm_,dat_,emb_,th_]:=Module[{links,medlink,mapp,me,con,li,pcms},
links=DeleteCases[Position[mcm,_?(#>th&)],_?(#[[1]]==#[[2]]&)];
medlink=Select[DeleteCases[Flatten[Outer[If[#1[[2]]==#2[[1]],Join[#1,{#2[[2]]}]]&,links,links,1],1],Null],
(MemberQ[links,{#[[1]],#[[-1]]}]&)];
mapp=(li=#;
me=Select[medlink,(#[[1]]==li[[1]]&&#[[-1]]==li[[2]]&)];
con={If[Length[me]==0,li,Join[{me[[1,1]],me[[1,-1]]},me[[All,2]]]]};
Which[Length[con[[1]]]==2,
li->Extract[mcm,li],
Length[con[[1]]]==3,
li->Pcm[dat[[con[[1]]]],emb[[{1,3}]]],
True,
pcms=(Pcm[dat[[#]],emb[[{1,3}]]])&/@Extcon[con];
While[Length[pcms]!=1,
pcms=Pcc/@Partition[pcms,3]];
li->pcms[[1]]]
)&/@links;
Fold[ReplacePart[#1,#2]&,Table[0,{Length[mcm]},{Length[mcm]}],mapp]
]

Mcm2Pcmp[mcm_,dat_,emb_,th_]:=Module[{links,medlink,mapp,me,con,li,pcms},
links=DeleteCases[Position[mcm,_?(#>th&)],_?(#[[1]]==#[[2]]&)];
medlink=Select[DeleteCases[Flatten[Outer[If[#1[[2]]==#2[[1]],Join[#1,{#2[[2]]}]]&,links,links,1],1],Null],
(MemberQ[links,{#[[1]],#[[-1]]}]&)];
mapp=ParallelTable[
me=Select[medlink,(#[[1]]==li[[1]]&&#[[-1]]==li[[2]]&)];
con={If[Length[me]==0,li,Join[{me[[1,1]],me[[1,-1]]},me[[All,2]]]]};
Which[Length[con[[1]]]==2,
li->Extract[mcm,li],
Length[con[[1]]]==3,
li->Pcm[dat[[con[[1]]]],emb[[{1,3}]]],
True,
pcms=(Pcm[dat[[#]],emb[[{1,3}]]])&/@Extcon[con];
While[Length[pcms]!=1,
pcms=Pcc/@Partition[pcms,3]];
li->pcms[[1]]],
{li,links},DistributedContexts->Automatic];
Fold[ReplacePart[#1,#2]&,Table[0,{Length[mcm]},{Length[mcm]}],mapp]
]

PCMall[dat_,th_]:=Module[{datx,emb,mcm,pcm},
datx=Transpose[dat];
emb=FEmbQ/@datx;
mcm=Table[Table[Mcm[datx[[y]],datx[[x]],emb[[y]]],{x,1,Length[dat[[1]]]}],{y,1,Length[dat[[1]]]}];
pcm=Mcm2Pcm[mcm,datx,emb,th];
{IdentityMatrix[Length[dat[[1]]]]+(pcm/._?(#<th&)->0),IdentityMatrix[Length[dat[[1]]]]+pcm(*tentative*)}
]

PCMallpara[dat_,th_]:=Module[{datx,emb,mcm,pcm},
datx=Transpose[dat];
emb=ParallelMap[FEmbQ[#]&,datx];
mcm=Table[Table[Mcm[datx[[y]],datx[[x]],emb[[y]]],{x,1,Length[dat[[1]]]}],{y,1,Length[dat[[1]]]}];
pcm=Mcm2Pcmp[mcm,datx,emb,th];
{IdentityMatrix[Length[dat[[1]]]]+(pcm/._?(#<th&)->0),IdentityMatrix[Length[dat[[1]]]]+pcm(*tentative*)}
]
(*********************************************************************************************************)


(*Correlation network*)
Corf[dat_,th_:0.01]:=Module[{corp,corpc},corp=Table[SpearmanRankTest[dat[[All,x]],dat[[All,y]]],{x,1,Length[dat[[1]]]},{y,1,Length[dat[[1]]]}];
corpc=corp/.{_?(#<th&)->1,_?(#>=th&)->0};
{SpearmanRho[dat]*corpc,corp}]

(**********************************************************************************************************)
SetVariables[data_,targetpos_,delay_]:=( 
{Drop[data,-delay],
Drop[data[[All,targetpos]],delay]})


RCsetup[pp_:{0.95(*spectral radius of reccurent matrix*),
0.95(*forgetting rate*),
0.001(*regularize coefficient in RLS*),
8(*Number of RLS update per one time step*)},
rep_:5000(*number of ESN population, if not given*)]:=( 

Unprotect[rcparameters];
Unprotect[p1];
Unprotect[p2];
Unprotect[p3];
Unprotect[p4];
Unprotect[q1];
Unprotect[q2];
Unprotect[InitializeNN];
Unprotect[comp];
Unprotect[compd];
Unprotect[RCcore];
Unprotect[RCprd];
Unprotect[RCmat];
Unprotect[RCbaselineprop];
Unprotect[RCcoren];
Unprotect[RCprdn];
Unprotect[RCmatn];
Unprotect[RCbaselinepropn];

{p1,p2,p3,p4}=pp;
q1=32;(*number of recurrent nodes*)
q2=rep;(*number of ESN population*)

rcparameters={p1,p2,p3,p4,q1,q2};

(***)
InitializeNN[nin_]:=Module[{nou,scwin,win,conn,crit,ww0,rho0,
gamma,ww,ilambda,delta,nx},
nx=q1;
scwin=0.1;
win=Table[scwin*RandomReal[{-1,1}],{nx},{nin}];
conn=0.1;
crit=1;
While[crit!=0,ww0=Table[If[RandomReal[]<conn,RandomReal[{-1,1}],0],{nx},{nx}];
crit=Count[Total[ww0]*Total[Transpose[ww0]],_?(Abs[#]<10^-8&)]];
rho0=Max[Abs[Eigenvalues[ww0]]];
gamma=p1*(1/rho0);
ww=gamma*ww0;
(**)
ilambda=1/p2;
delta=p3;
{nin,nx,win,ww,ilambda,delta}
];

comp=Compile[{{input,_Real,2},{output,_Real,1},
{nin,_Integer},{nx,_Integer},{win,_Real,2},{ww,_Real,2},{ilambda,_Real},{delta,_Real}},
Module[{inx,oux,xi,yi,pn,l1,ouprd,vn,gn,wou,yprd},
wou=Table[0.,{1},{nx}];
xi=Table[0.,{nx}];
pn=N[IdentityMatrix[nx]/delta];
yprd={0.};
MapThread[(inx=#1;oux={#2};
l1=Tanh[win . inx+ww . xi];
yprd=wou . l1;
wou=Nest[(ouprd=# . l1;
vn=oux-ouprd;
gn=Transpose[{ilambda*(pn . xi)/(1+ilambda*(xi . pn) . xi)}];
pn=ilambda*(pn-(gn . {xi})*pn);
(*state update*)
xi=l1;
#+1*Transpose[gn . {vn}])&,
wou,p4];
yprd[[1]])&,
{input,output}
]]];

compd=Compile[{{input,_Real,2},{output,_Real,1},
{nin,_Integer},{nx,_Integer},{win,_Real,2},{ww,_Real,2},{ilambda,_Real},{delta,_Real}},
Module[{inx,oux,xi,yi,pn,l1,ouprd,vn,gn,wou,yprd},
wou=Table[0.,{1},{nx}];
xi=Table[0.,{nx}];
pn=N[IdentityMatrix[nx]/delta];
yprd={0.};
MapThread[(inx=#1*RandomInteger[];oux={#2};
l1=Tanh[win . inx+ww . xi];
yprd=wou . l1;
wou=Nest[(ouprd=# . l1;
vn=oux-ouprd;
gn=Transpose[{ilambda*(pn . xi)/(1+ilambda*(xi . pn) . xi)}];
pn=ilambda*(pn-(gn . {xi})*pn);
(*state update*)
xi=l1;
#+1*Transpose[gn . {vn}])&,
wou,p4];
yprd[[1]])&,
{input,output}
]]];

RCcore[data_,targetpos_]:=Module[{itmax,input,target,ni,nl,baseprop,baggb,dummyin,tgest,
iactive,links,in,bagg0,pprev,eprev,aij,pij,newIactives,newins,baggs,inx,tgx,tgxe,
pbags,ebags,ebest,ppb,pbest,newiact,ppba},
itmax=q2;(*<-for reducing comp. time*)
{input,target}=SetVariables[data,targetpos,1];
ni=Length[First[data]];
nl=Length[data];
baseprop=RCbaselineprop[target];
(**)
iactive=ReplacePart[Table[0,{ni}],targetpos->1];(*initialize active variables*)
links=Total[iactive];
in=input[[All,Flatten[Position[iactive,1]]]];
bagg0=ParallelTable[
RMSE[target,Reservoirepd[in,target,InitializeNN[Length[in[[1]]]]]],
{itmax},DistributedContexts->Automatic];
pprev=Prop[bagg0];
eprev=Evp[pprev,baseprop];
aij=eprev[[1]]*iactive;
pij=(1-eprev[[2]])*iactive;
(**)
While[links<ni,
links=Total[iactive]+1;
newIactives=ReplacePart[iactive,#->1]&/@Flatten[Position[iactive,0]];
newins=input[[All,Flatten[Position[#,1]]]]&/@newIactives;
baggs=Transpose[ParallelTable[(in=#;
{inx,tgx}={in,target};
tgxe=Reservoirepd[inx,tgx,InitializeNN[Length[inx[[1]]]]];
RMSE[tgx,tgxe]
)&/@newins,
{itmax},DistributedContexts->Automatic]];
pbags=Prop/@baggs;
ebags=(Evp[#,pprev])&/@pbags;
ebest=Sort[ebags][[-1]];
ppb=Position[ebags[[All,1]],ebest[[1]]][[1,1]];
pbest=pbags[[ppb]];
newiact=newIactives[[ppb]];
ppba=Position[(newiact-iactive),1][[1,1]];
If[ebest[[1]]>0.&&ebags[[ppb]][[2]]<0.48,
iactive=newiact;
aij=ReplacePart[aij,ppba->ebest[[1]]];
pij=ReplacePart[pij,ppba->(1-ebest[[2]])];
pprev=pbest;
eprev=ebest,
Break[]
]];
{pprev,aij,1-pij,iactive}];

RCprd[data_,predyns_]:=Module[{iactives,itmax,ni,strps,targetpos,iactive,
input,target,in,bagg,tgtest},
iactives=predyns[[All,4]];
itmax=q2;
ni=Length[data[[1]]];
strps=Table[
targetpos=w;
iactive=iactives[[w]];
{input,target}=SetVariables[data,targetpos,1];
in=input[[All,Flatten[Position[iactive,1]]]];
bagg=ParallelTable[
tgtest=Reservoirep[in,target,InitializeNN[Length[in[[1]]]]];
{RMSE[target,tgtest],tgtest},{itmax},DistributedContexts->Automatic];
{Prop[bagg[[All,1]]],
(Quantile[#,{0.025,0.5,0.975}])&/@Transpose[bagg[[All,2]]],
target},
{w,1,ni}]];

RCmat[data_,predyns_]:=Module[{props,iactives,itmax,ni,strps,targetpos,prop,iactive,input,target,
iacs,iacz,newins,baggs,in,inx,tgx,tgxe,bagprop,evp},
props=predyns[[All,1]];
iactives=predyns[[All,4]];
itmax=q2;
ni=Length[data[[1]]];
strps=Table[
targetpos=w;
prop=props[[w]];
iactive=iactives[[w]];
{input,target}=SetVariables[data,targetpos,1];
iacs=(ReplacePart[iactive,#[[1]]->0]&)/@Position[iactive,1];
iacz=(ReplacePart[Table[0,{ni}],#[[1]]->1]&)/@Position[iactive,1];
newins=((#/.{}->{RandomReal[{0.8,1.2}]})&/@input[[All,Flatten[Position[#,1]]]])&/@iacs;
baggs=Transpose[ParallelTable[
(in=#;
{inx,tgx}={in,target};
tgxe=Reservoirepd[inx,tgx,InitializeNN[Length[in[[1]]]]];
RMSE[tgx,tgxe](*<- assuming a multivariate target//maxRMSE*)
)&/@newins,
{itmax},DistributedContexts->Automatic]];
bagprop=Prop/@baggs;
evp=Evp[#,prop]&/@bagprop;
{-Total[evp[[All,1]]* iacz],1-Total[evp[[All,2]]*iacz]},
{w,1,ni}];
Transpose[strps]];

RCbaselineprop[target_]:=Module[{dummyin,itmax,bagg,tgest},
itmax=q2;
bagg=ParallelTable[
dummyin=Table[{RandomReal[{0.8,1.2}]},{Length[target]}];
tgest=Reservoirepd[dummyin,target,InitializeNN[1]];
RMSE[target,tgest],
{itmax},DistributedContexts->Automatic];
Prop[bagg]];

(*non parallel*)

RCcoren[data_,targetpos_]:=Module[{itmax,input,target,ni,nl,baseprop,baggb,dummyin,tgest,
iactive,links,in,bagg0,pprev,eprev,aij,pij,newIactives,newins,baggs,inx,tgx,tgxe,
pbags,ebags,ebest,ppb,pbest,newiact,ppba},
itmax=q2;(*<-for reducing comp. time*)
{input,target}=SetVariables[data,targetpos,1];
ni=Length[First[data]];
nl=Length[data];
baseprop=RCbaselinepropn[target];
(**)
iactive=ReplacePart[Table[0,{ni}],targetpos->1];(*initialize active variables*)
links=Total[iactive];
in=input[[All,Flatten[Position[iactive,1]]]];
bagg0=Table[
RMSE[target,Reservoirepd[in,target,InitializeNN[Length[in[[1]]]]]],
{itmax}];
pprev=Prop[bagg0];
eprev=Evp[pprev,baseprop];
aij=eprev[[1]]*iactive;
pij=(1-eprev[[2]])*iactive;
(**)
While[links<ni,
links=Total[iactive]+1;
newIactives=ReplacePart[iactive,#->1]&/@Flatten[Position[iactive,0]];
newins=input[[All,Flatten[Position[#,1]]]]&/@newIactives;
baggs=Transpose[Table[(in=#;
{inx,tgx}={in,target};
tgxe=Reservoirepd[inx,tgx,InitializeNN[Length[inx[[1]]]]];
RMSE[tgx,tgxe]
)&/@newins,
{itmax}]];
pbags=Prop/@baggs;
ebags=(Evp[#,pprev])&/@pbags;
ebest=Sort[ebags][[-1]];
ppb=Position[ebags[[All,1]],ebest[[1]]][[1,1]];
pbest=pbags[[ppb]];
newiact=newIactives[[ppb]];
ppba=Position[(newiact-iactive),1][[1,1]];
If[ebest[[1]]>0.&&ebags[[ppb]][[2]]<0.48,
iactive=newiact;
aij=ReplacePart[aij,ppba->ebest[[1]]];
pij=ReplacePart[pij,ppba->(1-ebest[[2]])];
pprev=pbest;
eprev=ebest,
Break[]
]];
{pprev,aij,1-pij,iactive}];

RCprdn[data_,predyns_]:=Module[{iactives,itmax,ni,strps,targetpos,iactive,
input,target,in,bagg,tgtest},
iactives=predyns[[All,4]];
itmax=q2;
ni=Length[data[[1]]];
strps=Table[
targetpos=w;
iactive=iactives[[w]];
{input,target}=SetVariables[data,targetpos,1];
in=input[[All,Flatten[Position[iactive,1]]]];
bagg=Table[
tgtest=Reservoirep[in,target,InitializeNN[Length[in[[1]]]]];
{RMSE[target,tgtest],tgtest},{itmax}];
{Prop[bagg[[All,1]]],
(Quantile[#,{0.025,0.5,0.975}])&/@Transpose[bagg[[All,2]]],
target},
{w,1,ni}]];

RCmatn[data_,predyns_]:=Module[{props,iactives,itmax,ni,strps,targetpos,prop,iactive,input,target,
iacs,iacz,newins,baggs,in,inx,tgx,tgxe,bagprop,evp},
props=predyns[[All,1]];
iactives=predyns[[All,4]];
itmax=q2;
ni=Length[data[[1]]];
strps=Table[
targetpos=w;
prop=props[[w]];
iactive=iactives[[w]];
{input,target}=SetVariables[data,targetpos,1];
iacs=(ReplacePart[iactive,#[[1]]->0]&)/@Position[iactive,1];
iacz=(ReplacePart[Table[0,{ni}],#[[1]]->1]&)/@Position[iactive,1];
newins=((#/.{}->{RandomReal[{0.8,1.2}]})&/@input[[All,Flatten[Position[#,1]]]])&/@iacs;
baggs=Transpose[Table[
(in=#;
{inx,tgx}={in,target};
tgxe=Reservoirepd[inx,tgx,InitializeNN[Length[in[[1]]]]];
RMSE[tgx,tgxe](*<- assuming a multivariate target//maxRMSE*)
)&/@newins,
{itmax}]];
bagprop=Prop/@baggs;
evp=Evp[#,prop]&/@bagprop;
{-Total[evp[[All,1]]* iacz],1-Total[evp[[All,2]]*iacz]},
{w,1,ni}];
Transpose[strps]];

RCbaselinepropn[target_]:=Module[{dummyin,itmax,bagg,tgest},
itmax=q2;
bagg=Table[
dummyin=Table[{RandomReal[{0.8,1.2}]},{Length[target]}];
tgest=Reservoirepd[dummyin,target,InitializeNN[1]];
RMSE[target,tgest],
{itmax}];
Prop[bagg]];

Protect[rcparameters];
Protect[p1];
Protect[p2];
Protect[p3];
Protect[p4];
Protect[q1];
Protect[q2];
Protect[InitializeNN];
Protect[comp];
Protect[compd];
Protect[RCcore];
Protect[RCprd];
Protect[RCmat];
Protect[RCbaselineprop];
Protect[RCcoren];
Protect[RCprdn];
Protect[RCmatn];
Protect[RCbaselinepropn];
)

(************************************************************************************************************************************)

ES[p_,q_]:=(p[[1]]-q[[1]])/((((p[[3]]-1)*p[[2]]+(q[[3]]-1)*q[[2]])/(p[[3]]+q[[3]]-2))^0.5)

Prop[dist_]:={Median[dist],SmoothKernelDistribution[dist,"Silverman","Gaussian"]}

Reservoirep[inn_,out_,parameters_]:=Module[{nin,nx,win,ww,ilambda,delta,woui},
{nin,nx,win,ww,ilambda,delta}=parameters;
comp[inn,out,nin,nx,win,ww,ilambda,delta]
];

Reservoirepd[inn_,out_,parameters_]:=Module[{nin,nx,win,ww,ilambda,woui,delta},
{nin,nx,win,ww,ilambda,delta}=parameters;
compd[inn,out,nin,nx,win,ww,ilambda,delta]
];

BRMSE[target_,x_]:=Module[{cr,mtg,len},
(*,,\:30c8\:30e9\:30f3\:30b8\:30a7\:30f3\:30c8\:3092\:843d\:3068\:3059::10%??*)
cr=N[Count[target,Min[target]]/Length[target]];
len=Length[target];
If[cr>0.1,
mtg=Min[target];
Exp[-RootMeanSquare[If[#[[1]]<=mtg&&#[[2]]<#[[1]],0,#[[2]]-#[[1]]]&/@Take[Transpose[{target,x}],-IntegerPart[0.8*len]]]],
Exp[-RootMeanSquare[Take[x-target,-IntegerPart[0.8*len]]]]]
]

RMSE[ta_,tb_]:=Quiet[Exp[(-RootMeanSquare[Take[tb-ta,-IntegerPart[0.8*Length[tb]]]])]]

Evp[pa_,pb_]:={pa[[1]]-pb[[1]],1-CDF[pb[[2]]][pa[[1]]]}


(************************************************************************************************************************************)

(*misc*)

rndROCx[im_]:=Median[Table[ROCx[im,Partition[RandomSample[Flatten[im]],Length[im]],{0,1},Black],{16}]]

std[xx_]:=N[Transpose[Standardize[#]&/@Transpose[xx]]]
Norstd[xx_]:=N[Transpose[Standardize[#]&/@Transpose[xx]]]
Nordiffstd[xx_]:=N[Transpose[Standardize[#]&/@((Drop[RotateLeft[#]-#,-1])&/@Transpose[xx])]]
Logstd[xx_]:=N[Transpose[Standardize[#]&/@Transpose[Log[xx]]]]
Rnkstd[xx_]:=N[Transpose[(Standardize[Sort[MapIndexed[(Drop[Join[#1,#2],1])&,Sort[Transpose[{#,Range[Length[#]]}]]]][[All,2]]])&/@Transpose[xx]]]
LogRnkstd[xx_]:=N[Transpose[(Standardize[Sort[MapIndexed[(Drop[Join[#1,#2],1])&,Sort[Transpose[{#,Range[Length[#]]}]]]][[All,2]]])&/@Transpose[Log[xx]]]]

Rnkmat[mat_]:=Module[{rk,rma},rma=(mat)/.((#[[1]]->If[#[[1]]==0.,0,#[[2]]])&/@Transpose[{rk=Union[Flatten[mat]],Range[Length[rk]]}]);
N[rma/Max[rma]]]


Symm[mm_]:=MapThread[(Max[Abs[{#1,#2}]])&,{mm,Transpose[mm]},2]

Nev[act0_,est0_,rmode_,smode_,fmode_,linecol_:Black]:=Module[{procest,procact,uu,act,lrange,ppc,rocpoints,est,matched,precision},
If[smode==0,procest=Abs[est0];procact=Abs[act0],
procest=Symm[est0];
procact=Symm[act0]];
uu=DeleteCases[Union[Abs[Flatten[procest]]],Infinity];
act=procact/.{_?(Abs[#]<=10^-10&)->0,_?(Abs[#]>10^-10&)->1};
lrange=uu;
ppc=Sort[Table[If[rmode==0,
est=procest/.{_?(Abs[#]<=level&)->0,_?(Abs[#]>level&)->1},
est=procest/.{_?(Abs[#]<=level&)->1,_?(Abs[#]>level&)->0}];
matched=Transpose[(Flatten[MapIndexed[(Drop[#1,#2])&,#]])&/@{act,est}];
{1-N[Count[matched,_?(#[[1]]==0&&#[[2]]==0&)]/(Count[matched,_?(#[[1]]==0&)]/.(0->1))],
N[Count[matched,_?(#[[1]]==1&&#[[2]]==1&)]/(Count[matched,_?(#[[1]]==1&)]/.(0->1))],
N[Count[matched,_?((#[[1]]==1&&#[[2]]==1)||(#[[1]]==0&&#[[2]]==0)&)]/Length[matched]],
N[Count[matched,_?(#[[1]]==1&&#[[2]]==1&)]/(Count[matched,_?(#[[2]]==1&)]/.(0->1))],
N[Total[matched[[All,2]]]/Length[Flatten[matched]]]},
{level,If[rmode==0,lrange-10^-8,lrange+10^-8]}]];
rocpoints=Join[{{0,0}},Sort[ppc[[All,{1,2}]]],{{1,1}}];
{{Total[((#[[2,1]]-#[[1,1]])*#[[1,2]]+((#[[2,1]]-#[[1,1]])*(#[[2,2]]-#[[1,2]]))/2)&/@Partition[rocpoints,2,1]],(*ROCAUC*)
ppc[[All,3]][[-2]],(*Accuracy*)
(ppc[[All,4]])[[-2]]*(ppc[[All,2]])[[-2]]*2/((ppc[[All,2]])[[-2]]+(ppc[[All,4]]+10^-8)[[-2]]),(*f1-score*)
(ppc[[All,4]])[[-2]],(*precision*)
(ppc[[All,2]])[[-2]],(*sensitivity*)
1-(ppc[[All,1]])[[-2]](*specificity*)},
If[fmode==1,
rocpoints=Join[{{0,0}},Sort[ppc[[All,{1,2}]]],{{1,1}}];
Show[ListPlot[rocpoints,Joined->True,Mesh->All,AspectRatio->1,Joined->False,PlotRange->All,
MeshStyle->linecol,PlotStyle->linecol],
Plot[z,{z,0,1},PlotStyle->{Black,Dashed}],PlotRange->{{0,1},{0,1}},GridLines->Automatic,
Frame->True,FrameStyle->Directive[16,Black,FontFamily->"Helvetica"],FrameLabel->{"1-Specificity","Sensitivity"},ImageSize->400,Axes->False]]}]



End[]

(*LIMITS*)
Protect[LinearRegression];
Protect[GenerateXY];
Protect[GenerateXYn];
Protect[GenerateDdbDtg];
Protect[LIMITS];
Protect[LIMall];
Protect[LIMallpara];
(*CCM*)
Protect[Onestepah];
Protect[EmbQ];
Protect[Ccm];
Protect[Ccmrsmp];
Protect[CCMkendall];
Protect[CCMks];
Protect[CCMebisu];
Protect[CCMseason];
Protect[CCMallparasrg];
Protect[CCMallparanonsrg];
Protect[Seasonf]
Protect[DumFT];
Protect[EBZ];
Protect[DumSeason];
(*PCM*)
Protect[Falseneighbour];
Protect[FEmbQ];
Protect[Mcm];
Protect[Pcm];
Protect[Pcc];
Protect[SimplexProjection];
Protect[Extcon];
Protect[Mcm2Pcm];
Protect[Mcm2Pcmp];
Protect[PCMall];
Protect[PCMallpara];
(*Correkation*)
Protect[Corf];
(*EcohNet*)
Protect[SetVariables];
Protect[RCsetup];
Protect[ES];
Protect[Prop];
Protect[Reservoirep];
Protect[Reservoirepd];
Protect[BRMSE];
Protect[RMSE];
Protect[Evp];
Protect[rcparameters];
Protect[p1];
Protect[p2];
Protect[p3];
Protect[p4];
Protect[q1];
Protect[q2];
Protect[InitializeNN];
Protect[comp];
Protect[compd];
Protect[RCcore];
Protect[RCprd];
Protect[RCmat];
Protect[RCbaselineprop];
Protect[RCcoren];
Protect[RCprdn];
Protect[RCmatn];
Protect[RCbaselinepropn];
(*misc*)
Protect[rndROCx];
Unprotect[std];
Protect[Norstd];
Protect[Nordiffstd];
Protect[Logstd];
Protect[Rnkstd];
Protect[LogRnkstd];
Protect[Rnkmat];
Protect[Symm];
Protect[Nev];

EndPackage[]



