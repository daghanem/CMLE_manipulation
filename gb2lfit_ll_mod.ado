*! Shu Shen August 2016.
*! The program is modified from version 3.0.0 of gb2lfit by Stephen P. Jenkins


program define gb2lfit_ll_mod

        version 14.0
        args lnf ln_a ln_b ln_p ln_q
	
        quietly {

                tempname a b p q 
                scalar `a' = exp(`ln_a')
                scalar `b' = exp(`ln_b')
                scalar `p' = exp(`ln_p')
                scalar `q' = exp(`ln_q')
					if "$S_mlcens" == "" {
                        replace `lnf' = ln(`a') + (`a'*`p'-1)*ln($S_mlinc) ///
                                - `a'*`p'*ln(`b')   ///
                                - lngamma(`p') - lngamma(`q') + lngamma(`p'+`q') ///
                                - (`p'+`q')*ln(1+($S_mlinc/`b')^`a') 
					}
					
					if "$S_mlcens" != "" {
						tempvar lnd lnS 
                        ge double `lnd' = ln(`a') + (`a'*`p'-1)*ln($S_mlinc) ///
                                - `a'*`p'*ln(`b')   ///
                                - lngamma(`p') - lngamma(`q') + lngamma(`p'+`q') ///
                                - (`p'+`q')*ln(1+($S_mlinc/`b')^`a') 
						ge double `lnS'=ln( ibeta(`p',`q', (rtcff1/`b')^`a'/(1+(rtcff1/`b')^`a') ) - ibeta(`p',`q', ($S_mlinc/`b')^`a'/(1+($S_mlinc/`b')^`a') ) )							
						replace `lnf' = cond($S_mlcens, `lnS', `lnd', .)	
						
						if "$S_mlcens2" != "" {
							tempvar lnS2
							ge double `lnS2'=ln( ibeta(`p',`q', (rtcff2/`b')^`a'/(1+(rtcff2/`b')^`a') ) - ibeta(`p',`q', ($S_mlinc/`b')^`a'/(1+($S_mlinc/`b')^`a') ) )							
							replace `lnf' = cond($S_mlcens2, `lnS', `lnf', .)					 
						}
					}	
				}
		
		
		end
