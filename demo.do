clear all
cd "C:\Users\dghanem\Dropbox\Turning a Blind Eye\final\codes\CMLE_codes\DEMO"
set more off


capture program drop btstrp
program define btstrp, rclass
	tempvar Fxpeval Fxnpeval propFxbt polldayhatbt mu1bt mu2bt nu1bt Fxprteval nu2bt
	gb2lfit_mod pmconc, censvar(cens) from(intv)
	local esta = e(ba)
	local estb = e(bb)
	local estp = e(bp)
	local estq = e(bq)

	local atx=cff1
    gen `Fxpeval'=ibeta(`estp',`estq', (`atx'/`estb')^`esta'/(1+(`atx'/`estb')^`esta'))
	sum pmconc if pmconc<=float(cff1)
	gen `Fxnpeval'=r(N)/_N
	gen `propFxbt'=`Fxnpeval'-`Fxpeval'  // total proportion of manipulation up to the first cutoff
	sum `propFxbt'
	return scalar propFxbtstrp=r(mean)
	
	gen `polldayhatbt'=(1-`Fxpeval')*totday
	sum `polldayhatbt'
	return scalar polldayhatbtstrp=r(mean)

	gen `mu1bt'=`propFxbt'/`Fxnpeval'
	sum `mu1bt'
	return scalar mu1btstrp=r(mean)
	
	sum pmconc if pmconc<=float(lftcff1)
	gen `mu2bt'=`propFxbt'/(`Fxnpeval'-r(N)/_N)
	sum `mu2bt'
	return scalar mu2btstrp=r(mean)

	gen `nu1bt'=`propFxbt'/(1-`Fxpeval')
	sum `nu1bt'
	return scalar nu1btstrp=r(mean)

	local atx=rtcff1
    gen `Fxprteval'=ibeta(`estp',`estq', (`atx'/`estb')^`esta'/(1+(`atx'/`estb')^`esta'))
	gen `nu2bt'=`propFxbt'/(`Fxprteval'-`Fxpeval')
	sum `nu2bt'
	return scalar nu2btstrp=r(mean)
	
	return scalar estabtstrp=`esta'
	return scalar estbbtstrp=`estb'
	return scalar estpbtstrp=`estp'
	return scalar estqbtstrp=`estq'
	end
	
use "demodata.dta", clear

local idval=7 
local yearval=2006


scalar cff1=0.15
qui sum pmconc if pmconc<=float(cff1)
qui gen Fxnp=r(N)/_N  
qui gen cens=0


* manipulation window W1 
scalar lftcff1=0.135   
scalar rtcff1=0.22
replace cens=1 if pmconc>float(lftcff1) & pmconc <=float(cff1)  
gb2lfit_mod pmconc, cdf(pmparacdf1) censvar(cens)
local ahat=e(ba)
local bhat=e(bb)
local phat=e(bp)
local qhat=e(bq)
local atx=cff1
qui gen Fxp_w1=ibeta(`phat',`qhat', (`atx'/`bhat')^`ahat'/(1+(`atx'/`bhat')^`ahat'))
qui gen prop_w1=(Fxnp-Fxp_w1)/Fxnp
qui gen polldayhat_w1=(1-Fxp_w1)*totday
qui gen mu1_w1=prop_w1/Fxnp
qui sum pmconc if pmconc<=float(lftcff1)
qui gen mu2_w1=prop_w1/(Fxnp-r(N)/_N)
qui gen nu1_w1=prop_w1/(1-Fxp_w1)
local atx=rtcff1
qui gen Fxprt_w1=ibeta(`phat',`qhat', (`atx'/`bhat')^`ahat'/(1+(`atx'/`bhat')^`ahat'))
qui gen nu2_w1=prop_w1/(Fxprt_w1-Fxp_w1)

scalar lna0=ln(e(ba))
scalar lnb0=ln(e(bb))
scalar lnp0=ln(e(bp))
scalar lnq0=ln(e(bq))
matrix intv = (lna0, lnb0, lnp0, lnq0)
matrix colnames intv = ln_a:_cons ln_b:_cons ln_p:_cons ln_q:_cons

bootstrap ahatbs=r(estabtstrp) bhatbs=r(estbbtstrp) phatbs=r(estpbtstrp) qhatbs=r(estqbtstrp) ///
propbs=r(propFxbtstrp) polldayhatbs=r(polldayhatbtstrp) mu1bs=r(mu1btstrp)  mu2bs=r(mu2btstrp)  ///
nu1bs=r(nu1btstrp) nu2bs=r(nu2btstrp), cluster(woy) ///
reps(500) saving(bsreg,replace): btstrp 
append using bsreg.dta
local atx=cff1
qui egen prop_w1_bsse=sd(propbs)
qui egen polldayhat_w1_bsse=sd(polldayhatbs)
qui egen mu1_w1_bsse=sd(mu1bs)
qui egen mu2_w1_bsse=sd(mu2bs)
qui egen nu1_w1_bsse=sd(nu1bs)
qui egen nu2_w1_bsse=sd(nu2bs)

local cityname=city[1] 
log using log_demo, replace // documenting estimation results in Stata log file
* window 1, [0.135,0.22]
sum prop_w1 mu1_w1 mu2_w1 nu1_w1 nu2_w1 prop_w1_bsse mu1_w1_bsse mu2_w1_bsse nu1_w1_bsse nu2_w1_bsse 
log off

*manipulation window W2 
scalar lftcff1=0.135   
scalar rtcff1=0.18
replace cens=1 if pmconc>float(lftcff1) & pmconc <=float(cff1)  
gb2lfit_mod pmconc, cdf(pmparacdf2)  censvar(cens)
local ahat=e(ba)
local bhat=e(bb)
local phat=e(bp)
local qhat=e(bq)
local atx=cff1
qui gen Fxp_w2=ibeta(`phat',`qhat', (`atx'/`bhat')^`ahat'/(1+(`atx'/`bhat')^`ahat'))
qui gen prop_w2=(Fxnp-Fxp_w2)/Fxnp
qui gen polldayhat_w2=(1-Fxp_w2)*totday
qui gen mu1_w2=prop_w2/Fxnp
qui sum pmconc if pmconc<=float(lftcff1)
qui gen mu2_w2=prop_w2/(Fxnp-r(N)/_N)
qui gen nu1_w2=prop_w2/(1-Fxp_w2)
local atx=rtcff1
qui gen Fxprt_w2=ibeta(`phat',`qhat', (`atx'/`bhat')^`ahat'/(1+(`atx'/`bhat')^`ahat'))
qui gen nu2_w2=prop_w2/(Fxprt_w2-Fxp_w2)

scalar lna0=ln(e(ba))
scalar lnb0=ln(e(bb))
scalar lnp0=ln(e(bp))
scalar lnq0=ln(e(bq))
matrix intv = (lna0, lnb0, lnp0, lnq0)
matrix colnames intv = ln_a:_cons ln_b:_cons ln_p:_cons ln_q:_cons

bootstrap ahatbs=r(estabtstrp) bhatbs=r(estbbtstrp) phatbs=r(estpbtstrp) qhatbs=r(estqbtstrp) ///
propbs=r(propFxbtstrp) polldayhatbs=r(polldayhatbtstrp) mu1bs=r(mu1btstrp)  mu2bs=r(mu2btstrp)  ///
nu1bs=r(nu1btstrp) nu2bs=r(nu2btstrp), cluster(woy) ///
reps(500) saving(bsreg,replace): btstrp 
append using bsreg.dta
local atx=cff1
qui egen prop_w2_bsse=sd(propbs)
qui egen polldayhat_w2_bsse=sd(polldayhatbs)
qui egen mu1_w2_bsse=sd(mu1bs)
qui egen mu2_w2_bsse=sd(mu2bs)
qui egen nu1_w2_bsse=sd(nu1bs)
qui egen nu2_w2_bsse=sd(nu2bs)

local cityname=city[1] 
log on // documenting estimation results in Stata log file
* window 2, [0.135,0.18]
sum prop_w2 mu1_w2 mu2_w2 nu1_w2 nu2_w2 prop_w2_bsse mu1_w2_bsse mu2_w2_bsse nu1_w2_bsse nu2_w2_bsse 
log close


*****************
* POLYNOMIAL *
*****************
drop if pmconc==.
qui gen xrnd=round(pmconc,0.005) // round x to 0.005. Chetty's method requires discrete variables
gen xdisc=0.005*_n-0.005 
replace xdisc=. if xdisc>1
gen xfreq=.
forvalues i=1/201 {
sum xrnd if xrnd==float((`i'-1)*0.005)
replace xfreq=r(N)/_N in `i'
}
gen xdisc2=xdisc^2
gen xdisc3=xdisc^3
gen xdisc4=xdisc^4
gen xdisc5=xdisc^5
gen xdisc6=xdisc^6
gen xdisc7=xdisc^7

* order 6, window 1
forvalues i=28/45{
local val=(`i'-1)*5
gen xdum`val'=(xrnd==float((`i'-1)*0.005))
}
reg xfreq xdisc xdisc2 xdisc3 xdisc4 xdisc5 xdisc6 xdum*
forvalues i=28/45{
local val=(`i'-1)*5
replace xdum`val'=0
}
predict xfreqhat
sum xfreqhat in 1/31
gen Fxppoly61=r(mean)*r(N)
gen pmpolycdf61=sum(xfreqhat)
drop xdum* xfreqhat

* order 5, window 1
forvalues i=28/45{
local val=(`i'-1)*5
gen xdum`val'=(xrnd==float((`i'-1)*0.005))
}
reg xfreq xdisc xdisc2 xdisc3 xdisc4 xdisc5 xdum*
forvalues i=28/45{
local val=(`i'-1)*5
replace xdum`val'=0
}
predict xfreqhat
sum xfreqhat in 1/31
gen Fxppoly51=r(mean)*r(N)
gen pmpolycdf51=sum(xfreqhat)
drop xdum* xfreqhat

* order 7, window 1
forvalues i=28/45{
local val=(`i'-1)*5
gen xdum`val'=(xrnd==float((`i'-1)*0.005))
}
reg xfreq xdisc xdisc2 xdisc3 xdisc4 xdisc5 xdisc6 xdisc7 xdum*
forvalues i=28/45{
local val=(`i'-1)*5
replace xdum`val'=0
}
predict xfreqhat
sum xfreqhat in 1/31
gen Fxppoly71=r(mean)*r(N)
gen pmpolycdf71=sum(xfreqhat)
drop xdum* xfreqhat

* order 6, window 2
forvalues i=28/37{
local val=(`i'-1)*5
gen xdum`val'=(xrnd==float((`i'-1)*0.005))
}
reg xfreq xdisc xdisc2 xdisc3 xdisc4 xdisc5 xdisc6 xdum*
forvalues i=28/37{
local val=(`i'-1)*5
replace xdum`val'=0
}
predict xfreqhat
sum xfreqhat in 1/31
gen Fxppoly62=r(mean)*r(N)
gen pmpolycdf62=sum(xfreqhat)
drop xdum* xfreqhat

* order 5, window 2
forvalues i=28/37{
local val=(`i'-1)*5
gen xdum`val'=(xrnd==float((`i'-1)*0.005))
}
reg xfreq xdisc xdisc2 xdisc3 xdisc4 xdisc5 xdum135-xdum180
forvalues i=28/37{
local val=(`i'-1)*5
replace xdum`val'=0
}
predict xfreqhat
sum xfreqhat in 1/31
gen Fxppoly52=r(mean)*r(N)
gen pmpolycdf52=sum(xfreqhat)
drop xdum* xfreqhat

* order 7, window 2
forvalues i=28/37{
local val=(`i'-1)*5
gen xdum`val'=(xrnd==float((`i'-1)*0.005))
}
reg xfreq xdisc xdisc2 xdisc3 xdisc4 xdisc5 xdisc6 xdisc7 xdum*
forvalues i=28/37{
local val=(`i'-1)*5
replace xdum`val'=0
}
predict xfreqhat
sum xfreqhat in 1/31
gen Fxppoly72=r(mean)*r(N)
gen pmpolycdf72=sum(xfreqhat)
drop xdum* xfreqhat

cumul pmconc, gen(pmecdf) 
gen cutoff=pmconc*0+cff1

local cityname=city[1] 
line pmecdf pmconc if pmconc<=0.4, sort xlabel(0(0.1)0.4) lcolor(black) graphregion(color(white)) ///
|| line pmparacdf1 pmconc if pmconc<=0.4, lcolor(navy) sort ///
|| line pmparacdf2 pmconc if pmconc<=0.4, lcolor(navy) lpattern(dash) sort ///
|| line pmecdf cutoff, lcolor(black) lpattern(dot) title("`cityname', `yearval': Censored MLE")  ///
name(demo_MLE, replace) aspect(0.8) ///
legend(order(1 "ECDF" 2 "MLE-W1" 3 "MLE-W2") ///
rows(2) size(4) region(c(none)) bm(zero)) ytitle("") 

line pmecdf pmconc if pmconc<=0.4, sort xlabel(0(0.1)0.4) lcolor(black) graphregion(color(white)) ///
|| line pmpolycdf51 xdisc if xdisc<=0.4, lcolor(red) ///
|| line pmpolycdf61 xdisc if xdisc<=0.4, lcolor(orange) ///
|| line pmpolycdf71 xdisc if xdisc<=0.4, lcolor(purple) ///
|| line pmecdf cutoff, lcolor(black) lpattern(dot) title("`cityname', `yearval': Polynomial Fitting")  ///
name(demo_POLYw1, replace) aspect(0.8) ///
legend(order(1 "ECDF" 2 "Poly5-W1" 3 "Poly6-W1" 4 "Poly7-W1") ///
rows(2) size(4) region(c(none)) bm(zero)) ytitle("") 

line pmecdf pmconc if pmconc<=0.4, sort xlabel(0(0.1)0.4) lcolor(black) graphregion(color(white)) ///
|| line pmpolycdf52 xdisc if xdisc<=0.4, lcolor(red) lpattern(dash) ///
|| line pmpolycdf62 xdisc if xdisc<=0.4, lcolor(orange) lpattern(dash) ///
|| line pmpolycdf72 xdisc if xdisc<=0.4, lcolor(purple) lpattern(dash) ///
|| line pmecdf cutoff, lcolor(black) lpattern(dot) title("`cityname', `yearval': Polynomial Fitting")  ///
name(demo_POLYw2, replace) aspect(0.8) ///
legend(order(1 "ECDF" 2 "Poly5-W2" 3 "Poly6-W2" 4 "Poly7-W2" ) ///
rows(2) size(4) region(c(none)) bm(zero)) ytitle("") 

graph combine  demo_MLE  demo_POLYw1 demo_POLYw2, ///
rows(1) graphregion(color(white)) ysize(7) xsize(20)  
graph export "./outputgraphs/demo_compare.pdf", replace
