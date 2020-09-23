*! Shu Shen August 2016.
*! The program is modified from version 3.0.0 of gb2lfit by Stephen P. Jenkins


program define gb2lfit_mod, eclass 

version 14.0
        syntax varlist(max=1) [if] [in] [aw fw pw iw] [,  ///
                From(string) CENSvar(varname) CENSvar2(varname) ///
				CDF(namelist max=1) PDF(namelist max=1)  ///
                Robust Cluster(varname)  Level(integer $S_level) ///
                noLOG  * ]

        local title "Censored ML fit of GB2 distribution"

        local inc "`varlist'"
        

        if "`cdf'" ~= "" {
                confirm new variable `cdf' 
        }
        if "`pdf'" ~= "" {
                confirm new variable `pdf' 
        }



        local option0 `options'

        local wtype `weight'
        local wtexp `"`exp'"'
        if "`weight'" != "" { 
                local wgt `"[`weight'`exp']"'  
        }


        if "`weight'" == "pweight" | "`cluster'" != "" {
                local robust "robust"
        }

        if "`cluster'" ! = "" { 
                local clopt "cluster(`cluster')" 
        }

        if "`level'" != "" {
                local level "level(`level')"
        }
        
        local log = cond("`log'" == "", "noisily", "quietly") 

        marksample touse 
        markout `touse' `varlist' `cluster' `censvar' 
        mlopts mlopts, `options'

		if "`censvar'" != "" {
                capture quietly assert (`censvar' == 1 | `censvar' == 0) if `touse'
                if _rc != 0 {
                        di as error "censvar must be 0/1 variable"
                        exit 198
                }
        }

		if "`censvar2'" != "" {
                capture quietly assert (`censvar2' == 1 | `censvar2' == 0) if `touse'
                if _rc != 0 {
                        di as error "censvar2 must be 0/1 variable"
                        exit 198
                }
        }
		
        set more off

quietly {

        count if `inc' < 0 & `touse'
        local ct =  r(N) 
        if `ct' > 0 {
                noi di " "
                noi di as text "Warning: {res:`inc'} has `ct' values < 0;" _c
                noi di as text " not used in calculations"
          }

        count if `inc' == 0 & `touse'
        local ct =  r(N) 
        if `ct' > 0 {
                noi di " "
                noi di as text "Warning: {res:`inc'} has `ct' values = 0;" _c
                noi di as text " not used in calculations"
        }

        replace `touse' = 0 if `inc' <= 0


        count if `touse' 
        if r(N) == 0 {
                error 2000 
        }



        if "`from'" != ""  {
                local b0 "`from'"
        }

        global S_mlinc "`inc'"

		global S_mlcens "`censvar'"
		
		global S_mlcens2 "`censvar2'"

        `log' ml model lf gb2lfit_ll_mod (ln_a: ) (ln_b: ) (ln_p: ) (ln_q: )        ///
                `wgt' if `touse' , maximize  nocov iterate(1000) ///
                collinear title(`title') `robust' init(`b0')             ///
                search(on) `clopt' `level' `mlopts' `stdopts' `modopts'

        eret local depvar "`inc'"
		

        tempname b ba bb bp bq
        mat `b' = e(b)
        mat `ba' = `b'[1,"ln_a:"] 
        mat `bb' = `b'[1,"ln_b:"]
        mat `bp' = `b'[1,"ln_p:"]
        mat `bq' = `b'[1,"ln_q:"]

        eret matrix b_lna = `ba'
        eret matrix b_lnb = `bb'
        eret matrix b_lnp = `bp'
        eret matrix b_lnq = `bq'

                tempname e v            

                mat `e' = e(b)
                local a = exp(`e'[1,1])
                local b = exp(`e'[1,2])
                local p = exp(`e'[1,3])
                local q = exp(`e'[1,4])

                eret scalar ba = `a'
                eret scalar bb = `b'
                eret scalar bp = `p'
                eret scalar bq = `q'

                mat `v' = e(V)
                eret scalar se_a = `a' * sqrt(`v'[1,1])
                eret scalar se_b = `b' * sqrt(`v'[2,2])
                eret scalar se_p = `p' * sqrt(`v'[3,3])
                eret scalar se_q = `q' * sqrt(`v'[4,4])
       
                        /* Estimated GB2 c.d.f. */

                if "`cdf'" ~= "" {                      
                        qui ge `cdf' = ibeta(`p',`q', (`inc'/`b')^`a'/(1+(`inc'/`b')^`a') ) if `touse'
                        eret local cdfvar "`cdf'"
                }


                        /* Estimated GB2 p.d.f. */
        
                if "`pdf'" ~= "" {
                        qui ge `pdf' = (`a'*(`inc')^(`a'*`p'-1))*exp(lngamma(`p'+`q')) / (                ///
                                         (`b'^(`a'*`p'))*exp(lngamma(`p') + lngamma(`q'))  ///
                                          *( (1 +(`inc'/`b')^`a')^(`p'+`q') ) ///
                                        ) if `touse'
                        eret local pdfvar "`pdf'"
                }

		} // end -quietly block-

        Display, `level' `diopts'

end


program define Display

        syntax [,Level(int $S_level)  *]
        local diopts "`options'"
        if `level' < 10 | `level' > 99 {
                local level = 95
        }

        local plus "plus"


                ml display, level(`level') `diopts' `plus'
                _diparm ln_a, level(`level') exp prob label("Parameters a")
                di in smcl in gr "{hline 13}{c +}{hline 64}"
                _diparm ln_b, level(`level') exp prob label("b")
                di in smcl in gr "{hline 13}{c +}{hline 64}"
                _diparm ln_p, level(`level') exp prob label("p")
                di in smcl in gr "{hline 13}{c +}{hline 64}"
                _diparm ln_q, level(`level') exp prob label("q")
                di in smcl in gr "{hline 13}{c BT}{hline 64}"
		
		end


