import numpy as np
import matplotlib.pyplot as plt
def modalpoint(array:np.ndarray,eps=1e-6,verbose=False, _extra=0):
    # ;+
    # ;function modalpoint
    # ;	returns the mode of distribution defined by set of unbinned array values
    # ;
    # ;syntax
    # ;	arraymode=modalpoint(array,eps=eps,verbose=verbose)
    # ;
    # ;parameters
    # ;	array	[INPUT; required] array of values for which mode must be found
    # ;
    # ;keywords
    # ;	eps	[INPUT; default=1e-6] a small number
    # ;	verbose	[INPUT] controls chatter
    # ;	_extra	[JUNK] here only to prevent crashing the program
    # ;
    # ;description
    # ;	sort the array, divide into two sets around midpoint of range,
    # ;	choose the set that contains more elements, and repeat until
    # ;	the range is sufficiently small, and declare the midpoint of
    # ;	the range to be the mode
    # ;
    # ;warning
    # ;	may fail for multimodal distributions
    # ;
    # ;example
    # ;	for i=0,20 do print,i,modalpoint(randomn(seed,10000L)+i)
    # ;
    # ;history
    # ;	translated to IDL by Vinay Kashyap from C code written for BEHR
    # ;	  by Taeyoung Park c.2003 (MarMMVI)
    # ;	now correctly handles case when input is quantized (VK; SepMMVI)
    # ;	added edge case where if split is even, stops right there (VK; MayMMXXII)
    # ;-

# ;	usage
    ok='ok'
    # np=n_params()
    na=len(array)
    # if np eq 0 then ok='Insufficient parameters' else $
    if na == 0:
        ok='Input array is undefined'
    elif na < 2:
        ok='Array must have at least 2 elements'
    if ok != 'ok':
        print,'Usage: arraymode=modalpoint(array,eps=eps,verbose=verbose)'
        print,'  return mode of array'
    # if np ne 0 then message,ok,/informational
        return np.nan

    # ;	inputs and some special cases
    if na < 3:
        return np.mean(array)
    ok=np.argwhere(np.isfinite(array)).flatten()
    ok = np.array(ok,dtype=int)
    mok = np.sum(np.isfinite(array))
    if mok == 0:
        return np.nan
    # print(ok)
    arr=np.take(array,ok) 
    # ok = np.array(ok,dtype=np.int64)
    os=np.sort(arr) 
    arr=os
    # ;
    vv=0
    if verbose:
        vv=verbose
    # ;
    # if keyword_set(eps) then epsilon=double(eps[0]) else epsilon=1d-6

    # ;	step through the array and find the mode
    go_on=1
    narr=len(arr) 
    amax=np.nanmax(arr)
    amin=np.nanmin(arr)
    while go_on:
        # if vv > 10:
        #     print(strtrim(narr,2)+'.. ',format='($,a)')
        # o1=where(arr gt 0.5*(amin+amax),mo1)
        o1 = np.argwhere(arr>0.5*(amin+amax)).flatten()
        mo1=len(o1)
        if mo1 == 0 or mo1 == narr:
            # ;message,'BUG?'
            # ;not a bug, this means that they are all identical, so quit right here
            return np.nanmedian(arr)

        if o1[0] > narr/2:
            tmparr=arr[0:o1[0]] 
        else:
            tmparr=arr[o1]
        if o1[0] == narr/2:
            # if vv gt 0 then message,'evenly split, might be multimodal; cannot deal',/informational
            tmparr=arr
            go_on=0	#;quit at this stage

        # if vv gt 100 then print,narr/2,mo1,o1[0],min(tmparr),max(tmparr)
        arr=tmparr
        narr=len(arr)
        amax=np.nanmax(arr)
        amin=np.nanmin(arr)

        if narr == 1:
            go_on=0	#;stop when there is only one element
        if amax-amin < eps:
            go_on=0	#;stop when range gets too small
    # endwhile

    return 0.5*(amin+amax)

if __name__=='__main__':

    arr=[1,2,3,3,2,5,6,6,7,6,6,8,8,8,7,1,2,1,9]

    md=modalpoint(arr)
    print(md)

    plt.hist(arr)
    plt.show()