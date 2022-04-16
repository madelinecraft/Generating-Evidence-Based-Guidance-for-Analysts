libname s10_50 "/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data10_50/output";

*proc import datafile="/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data10_50/data10_50.txt" out=s10_50.data dbms = dlm replace; *run;

data S10_50.DATA;
	%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
	infile '/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data10_50/data10_50.txt' delimiter = ' '
	MISSOVER DSD lrecl=32767 firstobs=2 ;
	informat replication $4. ;
	informat _imp best32. ;
	informat _id best32. ;
	informat subj best32. ;
	informat W2_b best32. ;
	informat W1 best32. ;
	informat time $3. ;
	informat y best32. ;
	informat X2_b best32. ;
	informat X1 best32. ;
	format replication $4. ;
	format _imp best12. ;
	format _id best12. ;
	format subj best12. ;
	format W2_b best12. ;
	format W1 best12. ;
	format time $3. ;
	format y best12. ;
	format X2_b best12. ;
	format X1 best12. ;
	input
	replication $
	_imp
	_id
	subj
	W2_b
	W1
	time $
	y
	X2_b
	X1;
	if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
run;

* convert replication to numeric;
data s10_50.data;
	set s10_50.data;
    rep_n = input(replication, 4.);
run;

* sort data;
proc sort data = s10_50.data; by rep_n _imp subj time; run;

* rename _imp variable for proc mianalyze;
data s10_50.data; rename _imp=_imputation_; set s10_50.data; run;

* grand-mean centering W1; 
data s10_50.data; set s10_50.data; 
W1_c = W1; *copy the variable in a datastep so you don't overwrite it in the next step;
run; 
proc standard data=s10_50.data out=s10_50.data mean=0; *mean=0 does the centering; 
by rep_n _imputation_;
var W1_c; 
run; 

* person-mean centering X1; 
proc means data=s10_50.data nway noprint; 
by rep_n _imputation_;
class subj; 
var X1;
output out=s10_50.means1 /*this makes a new dataset with only 1 row per person (per rep & imp) that contains the mean of X1 for that person*/ 
mean=X1_mean;  
run; 
proc sort data = s10_50.data; by rep_n _imputation_ subj time; run;
proc sort data=s10_50.means1; by rep_n _imputation_ subj; run; 
data s10_50.data;
merge s10_50.data s10_50.means1;
by rep_n _imputation_ subj;  
X1_c = X1-X1_mean; *center on the person-means, cw stands for "centered within";
run; 

* override "SAS set option OBS=0 and will continue to check statements";
options obs=max;
options replace;
options nosyntaxcheck;
OPTIONS NODMSSYNCHK;

* fit the model by replication and imputation;
ods graphics off; * suppress ods graphics for a more efficient simulation;
ods exclude all; * suppress ods output for a more efficient simulation;
ods output ConvergenceStatus = s10_50.myconvstat
		   ParameterEstimates = s10_50.myparamest
		   CovMatParmEst = s10_50.mycovmat 
		   FitStatistics = s10_50.myfitstats
		   Dimensions = s10_50.mydimensions; * output ods tables of interest to datasets;
proc nlmixed data=s10_50.data gconv = 1e-12 cov;
parms b0 = -.0637 b1 = .1002 b2 = .1426 b3 = .1015 b4 = .1188 lp0 = .1936 tau0 = .3595
lp1 = .0543 lp2 = .1352 tau1 = .0979 tau2 = .1403 tau3 = .0847 tau4 = .0972 var2 = 0 vcov12 = 0;
BY rep_n _imputation_;
*WHERE rep_n < 11 and _imputation_ < 11;
*where rep_n = 5 and _imputation_ = 5;
vare = exp(tau0 + X1_c*tau1 + X2_b*tau2 + W1_c*tau3 + W2_b*tau4 + u2);
varu = exp(lp0 + W1_c*lp1 + W2_b*lp2);
z = b0 + b1*X1_c + b2*X2_b + b3*W1_c + b4*W2_b + u1;
model y ~ normal(z,vare);
random u1 u2 ~ normal([0,0], [varu,vcov12,var2]) subject=subj;
*random u1 ~ normal([0], [varu]) subject=subj;
run;
ods output close;
ods exclude none; 

* override "SAS set option OBS=0 and will continue to check statements";
options obs=max;
options replace;
options nosyntaxcheck;
OPTIONS NODMSSYNCHK;

* create an indicator of non-convergent solutions for myparamest;
data s10_50.myparamest;
	set s10_50.myparamest;
	ind = .;
	if Gradient > 1e-2 or Gradient < -1e-2 or StandardError = . or StandardError = 0 then ind = 1;
	else ind = 0;
run;

* if indicated, make Estimate missing;
data s10_50.myparamest2; set s10_50.myparamest;
if ind = 1 then Estimate = .;
run; 

* sort by replication;
proc sort data = s10_50.myparamest2; by rep_n; run;

* remove missing estimates;
data s10_50.myparamest3;
  do until (last.rep_n);
    set s10_50.myparamest2;
    by rep_n;
    anymissing=max(anymissing,cmiss(of Estimate));
  end;
  do until (last.rep_n);
    set s10_50.myparamest2;
    by rep_n;
    if not anymissing then output;
  end;
run;

* combine parameter estimates with covariance matrix;
proc sort data = s10_50.mycovmat; by rep_n _imputation_; run;
proc sort data = s10_50.myparamest3; by rep_n _imputation_; run;
data s10_50.alloutput;
merge s10_50.myparamest3 s10_50.mycovmat;
by rep_n _imputation_;
run;

* delete any row with missing estimates;
data s10_50.alloutput2; set s10_50.alloutput;
if Estimate = "." then delete;
run;

* separate paramest estimates and covariance matrix;
data s10_50.myparamest4; set s10_50.alloutput2 (keep = rep_n _imputation_ Parameter Estimate Gradient StandardError DF tValue Probt Alpha Lower Upper); run;
data s10_50.mycovmat2; set s10_50.alloutput2 (keep = rep_n _imputation_ Parameter Row b0 b1 b2 b3 b4 lp0 tau0 lp1 lp2 tau1 tau2 tau3 tau4 var2 vcov12); run;

* re-number replication and imputation;
data s10_50.numbering; set s10_50.myparamest (keep = rep_n _imputation_); run;
data s10_50.myparamest5; merge s10_50.myparamest4 s10_50.numbering; 
if Estimate = "." then delete; run;
data s10_50.mycovmat3; merge s10_50.mycovmat2 s10_50.numbering; 
if Row = "." then delete; run;

* pool estimates across imputations;
ods graphics off; * suppress ods graphics for a more efficient simulation;
ods exclude all; * suppress ods output for a more efficient simulation;
* pool estimates across imputations for each replication;
* proc mianalyze resource: file:///C:/Users/mcraft/Downloads/2005_Chapter_IncompleteDataAndSAS%20(1).pdf;
ods output	ParameterEstimates = s10_50.pooledparamest
			VarianceInfo = s10_50.pooledvarinfo; * output ods tables of interest to datasets;
proc mianalyze parms=s10_50.myparamest5 covb = s10_50.mycovmat3 wcov bcov tcov;
by rep_n;
modeleffects b0 b1 b2 b3 b4 lp0 tau0 lp1 lp2 tau1 tau2 tau3 tau4 var2 vcov12;
run;
ods output close; * close dataset output;
ods exclude none; * unsuppress output;

* create flag for pooled datasets with zero between-imputation variance;
data s10_50.pooledvarinfo2; set s10_50.pooledvarinfo; 
flag = 0;
if BetVar = 0 then flag = 1;
run;
data s10_50.pooledparamest2; set s10_50.pooledparamest; 
flag = 0;
if Probt = . then flag = 1;
run;

* remove flagged rows;
data s10_50.pooledvarinfo3; set s10_50.pooledvarinfo2;
if flag = 1 then delete; run;

data s10_50.pooledparamest3; set s10_50.pooledparamest2;
if flag = 1 then delete; run;

* summarize pooled estimates;
data s10_50.pooledparamest3; set s10_50.pooledparamest3; by rep_n; 

if parm = 'b0' then trueval    = 0;
if parm = 'b1' then trueval    = .2;
if parm = 'b2' then trueval    = .1;
if parm = 'b3' then trueval    = .1;
if parm = 'b4' then trueval    = .1;
if parm = 'tau0' then trueval = .2;
if parm = 'lp1' then trueval    = .2;
if parm = 'lp2' then trueval    = .1;
if parm = 'tau1' then trueval    = .1;
if parm = 'tau2' then trueval    = .2;
if parm = 'tau3' then trueval    = .2;
if parm = 'tau4' then trueval    = .1;

if parm = 'lp0' then trueval  = 0.2;
if parm = 'var2' then trueval  = 1.1;
if parm = 'vcov12' then trueval = 0.3;

inrange = 0;
if (trueval > LCLMean and trueval < UCLMean) then inrange = 1;

width = UCLMean - LCLMean;

proc sort; by parm;

proc means noprint; var estimate trueval inrange width;  
output out = s10_50.pooledparamest4  mean(estimate trueval width) = estm truev widthm
                                  sum(inrange) = numin std(estimate) = estsd
						          n(estimate) = nconv;
by parm;
run;

data s10_50.pooledparamest5; set s10_50.pooledparamest4;
coverage = numin / nconv;
bias = estm - truev;
stdbias  = 100 * (bias/estsd);
rmse     = sqrt((estm - truev)**2 + estsd*estsd);

proc print; var parm truev nconv estm bias stdbias coverage rmse widthm;
run;

/*/
* write out the results;
DATA L30_500.out2; set L30_500.resultls;
FILE 'C:/Users/mcraft/Desktop/Dissertation Simulation/p_0/sim0.2_1.1_0.3/sim30_500/results0.2_1.1_0.3.txt';
put parameter $ 11-19 (truev) (6.3) (nconv) (5.0) 
    (estm bias stdbias) (12.6) (coverage) (8.4) (rmse widthm) (12.6);
run;
/*/
