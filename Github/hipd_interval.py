fsample = None
clev = 0.68
pdfnorm = None
fmode = None
verbose = None
e=0

import numpy as np
from mode import modalpoint
def hipd_interval(f,x=None,fsample=fsample,clev=clev,pdfnorm=pdfnorm,fmode=fmode,verbose=verbose, _extra=e):
# ;function	hipd_interval
# ;	computes and returns the interval [lower_bound,upper_bound] at a
# ;	specified confidence level that includes the highest probability
# ;	densities.  by definition, this is the shortest possible interval.
# ;
# ;syntax
# ;	hpd=hipd_interval(f,x,/fsample,clev=clev,pdfnorm=pdfnorm,$
# ;	fmode=fmode,verbose=verbose)
# ;
# ;parameters
# ;	f	[INPUT; required] the array for which the confidence interval
# ;		must be computed
# ;		* assumed to be a density function unless FSAMPLE is set
# ;	x	[INPUT; optional] abscissae values NOT FUNCTIONAL
# ;		* if not given, and F is a density function, then taken
# ;		  to be the array indices
# ;		* ignored if FSAMPLE is set
# ;		* if N(X).GT.1 but N(X).NE.N(F), X is ignored and
# ;		  FSAMPLE is set automatically
# ;
# ;keywords
# ;	fsample	[INPUT] if set, assumes that F is a set of samples from a
# ;		density function, as opposed to being the density function
# ;		itself NOT FUNCTIONAL NEEDS TO BE FUNCTIONAL
# ;	clev	[INPUT] confidence level at which to compute the intervals
# ;		* default is 0.68
# ;		* if < 0, abs(CLEV) is used
# ;		* if > 1 and < 100, then assumed to be given as a percentage
# ;		* if > 100, then 1-1/CLEV is used
# ;	pdfnorm	[INPUT] if set, forces F to integrate to abs(PDFNORM)
# ;		* if explicitly set to 0, does not normalize F at all
# ;		* if not set, normalizes to 1
# ;		* ignored if FSAMPLE is set
# ;		* WARNING: do not use this keyword unless you know
# ;		  what you are doing
# ;	fmode	[OUTPUT] the mode of the distribution
# ;	verbose	[INPUT] controls chatter
# ;	_extra	[INPUT ONLY] pass defined keywords to subroutines
# ;		MODALPOINT: EPS
# ;
# ;subroutines
# ;	MODALPOINT
# ;
# ;description
# ;	* if density function, find the cumulative integral around the
# ;	  mode and interpolate at CLEV to derive the HPD interval
# ;	* if array of sample values, find all intervals corresponding to
# ;	  CLEV that contain the mode and pick the smallest of the lot.
# ;	  this is a method devised by Vinay K and Taeyoung Park during
# ;	  the course of developing BEHR.
# ;	* works well only for unimodal distributions (and for multimodal
# ;	  distributions where all the modes are within the interval),
# ;	  but that's better than nothing.
# ;	* note that this interval is not invariant under transformations.
# ;	  for that, one must use equal-tail intervals, see eqt_interval()
# ;
# ;example
# ;	for i=1,20 do print,hipd_interval(randomn(seed,10000L)*i,/fsample)
# ;
# ;history
# ;	vinay kashyap (Mar2006)
# ;	bug correction with FSAMPLE cdf (VK; Apr2006)
# ;	added keyword FMODE (VK; Nov2006)
# ;	bug correction with F(X) case (VK; Apr2007)
# ;	now handles NaNs in input (VK; Apr2014)
# ;-

# ;	usage
    ok='ok' 
    # np=n_params()
    nf=len(f)
    # print(nf)
    if x:
        nx=len(x)
    else:
        nx=None
    # if np eq 0 then ok='Insufficient parameters' else $
    #  if nf eq 0 then ok='F is not defined' else $
    #   if nf lt 2 then ok='array size too small'
    # if ok ne 'ok' then begin
    # print,'Usage: hpd=hipd_interval(f,x,/fsample,clev=clev,pdfnorm=pdfnorm,$'
    # print,'       fmode=fmode,verbose=verbose)'
    # print,'  compute highest probability density interval'
    # #   if np ne 0 then message,ok,/informational
    # #   return,-1L
    # endif

    # ;	figure out density function or array
    dens=1
    if (fsample is not None):
        dens=0
    if nx:
        if nx > 1 and nx != nf:
            dens=0
    # ;
    if (dens):
        xx=np.arange(nf) 
        # print(xx)
        if nx:
            if nx == nf:
                xx=x
    # ;
    ff=f

    # ;	keywords
    vv=0
    if (verbose):
        vv=np.long(verbose[0]) > 1
    # ;
    crlev=0.68
    if (clev):
        crlev=0.0+clev
    if crlev < 0:
        crlev=abs(crlev)
    if crlev >= 1 and crlev < 100:
        crlev=crlev/100.
    if crlev >= 100:
        crlev = 1.0 - 1.0/crlev
    # print(crlev)
    # ;	find the mode
    fmax=np.nanmax(ff)
    imx = np.nanargmax(ff)
    fmin = np.nanmin(ff)
    # print(ff)
    # print(fmax,fmin,imx)
    if (dens):
        fmode=xx[imx]
        # print(fmode)
    else:
        fmode=modalpoint(ff,verbose=vv, _extra=e)

    # ;	compute interval
    if (dens):		#;(if prob density
        # ;sort
        os=np.argsort(xx) 
        # print(os)
        xx=xx[os.astype(np.int32)] 
        ff = np.array(ff)
        ff=ff[os.astype(np.int32)]
        # print(ff)
        #;get cdf
        cff = np.cumsum(ff)
        ii=np.arange(nf,dtype=np.int32)
            #;normalize to 1, or not, or to PDFNORM
        if (pdfnorm) is None:
            cff=cff/np.nanmax(cff)
            print(cff)
        else:
            if (pdfnorm):
                cff=np.abs(pdfnorm[0])*cff/np.nanmax(cff,)
        # endelse
            # ;this is where the mode lies on the cdf
        cfmode=np.interp(fmode,xx,cff)
            # ;given the mode's level, how far back can we go?
        cfmin=cfmode-crlev 
        cfmin = cfmin if cfmin>0 else 0
            # ;and how far up can we go?
        cfmax=cfmode+crlev 
        cfmax = cfmax if cfmax<1 else 1

            # ;the cdf levels translate to these bin indices
        ixmode=np.int32(np.interp(cfmode,cff,ii)) 
        ixmode = ixmode if ixmode > 0 else 0 
        ixmin=np.int32(np.interp(cfmin,cff,ii))
        ixmin = ixmin if ixmin > 0 else 0 
        # print(ixmin,ixmode)
            # ;set up the loop to find the smallest range
        go_on=1 
        k=ixmin 
        cfmax=cff[k]+crlev 
        drng0=np.nanmax(xx)-np.nanmin(xx,)
        hpdm=[]
        hpdp=[]
        while go_on:		#;{check every interval and pick the smallest
            # ;the highest index given the current lower index
            ixmax=np.int32(np.interp(cfmin+crlev,cff,ii)) 
            # print(cfmin,crlev,cff,ii)
            ixmax = ixmax if ixmax < (nf-1) else nf-1
            # ;and the range that corresponds to these indices
            drng=np.abs(xx[k]-xx[ixmax])
            # ;check if current interval is smaller
            # print(k)
            # print(drng0)
            if drng < drng0:
                # ;update smallest interval
                # print('here')
                hpdm=xx[k] 
                hpdp=xx[ixmax] 
                drng0=drng
                if vv > 1000:
                    print(k,hpdm,hpdp,drng0,cff[k],cff[ixmax],cfmode)
            			#;DRNG<DRNG0)
            # ;next step
            k=k+1
            cfmax=cff[k]+crlev
            if cff[k] > cfmode:
                go_on=0	#;quit if hit the mode
            if cfmax >= 1:
                go_on=0		#;quit if bumped up to the end
            if vv > 500 and go_on == 0:
                pass
                #idk
                # stop,'HALTing; type .CON to continue'
        			#;GO_ON}
           # ;and done
        return [hpdm,hpdp]

        # ;	;need the reverse() because we want to start from the peak
        # ;os=reverse(sort(ff))
        # ;	;subtract fmin to account for cases where ff drops below zero
        # ;if fmin lt 0 then cff=total(ff[os]-fmin,/cumul) else $
        # ;	cff=total(ff[os],/cumul)
        # ;	;normalize to 1, or not, or to PDFNORM
        # ;if not arg_present(pdfnorm) then cff=cff/max(cff) else begin
        # ;  if keyword_set(pdfnorm) then cff=abs(pdfnorm[0])*cff/max(cff)
        # ;endelse
        # ;	;sort the indices
        # ;xx=xx[os]
        # ;	;pick out those indices that fall on either side of mode
        # ;om=where(xx le fmode,mom) & op=where(xx ge fmode,mop)
        # ;	;in case mode is at extreme, peg the range to mode
        # ;hpdm=fmode & hpdp=fmode
        # ;	;interpolate on integrated function
        # ;if mom gt 1 then hpdm=interpol(xx[om],cff[om],crlev)
        # ;if mop gt 1 then hpdp=interpol(xx[op],cff[op],crlev)
        # ;	;and done
        # ;if vv gt 100 then print,hpdm,hpdp,fmode
        # ;if vv gt 500 then stop,'HALTing; type .CON to continue'
        # ;return,[hpdm,hpdp]
    else:			#;DENS)(if array of values #NOT USED?
        # ;first make sure everything is sorted
        os=np.argsort(ff)
        ii=np.arange(nf,dtype=np.int32)
            # ;get cdf
            # ;no PDFNORM nonsense with samples
        cff=np.arange(nf,dtype=np.double)/(nf-1.)
            # ;the mode is at this cumulative level
        # print(ff,os)
        ff=np.array(ff)
        cfmode=np.interp(fmode,ff[os],cff)
            # ;given the mode's level, how far back can we go?
        
        cfmin=cfmode-crlev 
        cfmin = cfmin if cfmin>0 else 0
            # ;and how far up can we go?
        cfmax=cfmode+crlev 
        cfmax = cfmax if cfmax<1 else 1

            # ;the cdf levels translate to these bin indices
        ixmode=np.int32(np.interp(cfmode,cff,ii)) 
        ixmode = ixmode if ixmode > 0 else 0 
        ixmin=np.int32(np.interp(cfmin,cff,ii))
        ixmin = ixmin if ixmin > 0 else 0 
        
        go_on=1 
        k=ixmin 
        cfmax=cff[k]+crlev 
        drng0=fmax-fmin
        hpdm=np.nan
        hpdp=np.nan
        while go_on:		#;{check every interval and pick the smallest
            # ;the highest index given the current lower index
            ixmax=k+np.int32(crlev*nf+0.5)
            # print(cfmin,crlev,cff,ii)
            ixmax = ixmax if ixmax < (nf-1) else nf-1
            # ;and the range that corresponds to these indices
            drng=abs(ff[os[k]]-ff[os[ixmax]])
            # ;check if current interval is smaller
            if drng <= drng0:		#;(found smaller interval
                # ;update smallest interval
                hpdm=ff[os[k]] 
                hpdp=ff[os[ixmax]] 
                drng0=drng
                if vv > 100:
                    print(k,hpdm,hpdp,drng0,cff[k],cfmode)
            			#;DRNG<DRNG0)
            # ;next step
            k=k+1
            cfmax=cff[k]+crlev
            if cff[k] > cfmode:
                go_on=0	#;quit if hit the mode
            if cfmax >= 1:
                go_on=0		#;quit if bumped up to the end
            if vv > 500 and go_on == 0:
                pass
                #idk
                # stop,'HALTing; type .CON to continue'
        			#;GO_ON}
            # ;and done
        return [hpdm,hpdp]
    					#;not DENS)

    # return !values.F_NAN	;should never get here
# print(__name__)
import matplotlib.pyplot as plt
from arviz import hdi
if __name__ == '__main__':
    arr=[1,1,1,22,32,1,22,1,4,5,6,74,2,1,1,1,1,2,2,2,2,2,2,2,2]
    hpd = hipd_interval(arr,clev=0.68,fsample=True)
    hdi_a = hdi(np.array(arr),hdi_prob=0.68)
    print(hpd,hdi_a)
    plt.hist(arr,bins=100)
    # plt.axvspan(arr[hpd[0]],arr[hpd[1]],alpha=0.2,color='purple')
    # plt.axvspan(hpd[0],hpd[1],alpha=0.2,color='purple')
    # plt.axvspan(hdi_a[0],hdi_a[1],alpha=0.2,color='red')
    plt.show()