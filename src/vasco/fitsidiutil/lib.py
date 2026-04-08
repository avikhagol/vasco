import numpy as np
from astropy.io import fits
from typing import List

class IdiHDU(fits.PrimaryHDU):
    @classmethod
    def match_header(cls, header):
        try:
            keyword = header.cards[0].keyword
        except:
            keyword = header.ascard[0].key
            pass
        
        return (keyword == 'SIMPLE' and 'GROUPS' in header and
            header['GROUPS'] == True and 'NAXIS' in header and header['NAXIS'] == 0)

fits.register_hdu(IdiHDU)

def dict_baseline(fitsfile=None,hdul=None):
    """
    for each baseline stores tuple of (distance,label,(id_an1,id_an2))
    """
    from scipy.spatial import distance
    if hdul is None:
        if fitsfile is None: raise ValueError("missing fits file path")
        else:hdul=fits.open(fitsfile)
    xyz,annames=hdul['ARRAY_GEOMETRY'].data.STABXYZ, hdul['ARRAY_GEOMETRY'].data.ANNAME

    dict_baseline={}
    
    for i,an1 in enumerate(annames):
        for j,an2 in enumerate(annames):
                d=distance.euclidean(xyz[i],xyz[j])*.001
                baseline_label=f"{an1}-{an2}"
                baseline_id=(i+1)*256+(j+1)
                dict_baseline[baseline_id]=(d,baseline_label,(i+1,j+1))
    return dict_baseline

def _gethduname(hdul, hdu_names:List[str]):
    hdu_id = []
    hdu_name = ''
    for i,hdu_found in enumerate(hdul.names) :
        for col_given in hdu_names:
            if col_given in hdu_found:
                hdu_id.append(i) 
                hdu_name = hdu_found
    return [hdu_name, hdu_id]

def _getcolname(hdu, cols:List[str]):
    return [col_found for col_found in hdu.cols for col_given in cols if col_given in col_found][0]


def get_yyyymmdd(dateobs):
    """
    Takes FITS DATE-OBS with fromat : yy/mm/dd or dd-mm-yyyy and returns (yyyy,mm,dd)
    Args:
        dateobs (str): FITS DATE-OBS/RDATE format

    Returns:
        tuple: (yyyy,mm,dd)
    """    
    yyyy, mm, dd = 0,0,0
    try:
        dateobs = [int(d) for d in dateobs.split('-')]
        yyyy = dateobs[0]
        mm, dd = dateobs[1], dateobs[2]
    except:
        dateobs = [int(d) for d in dateobs.split('/')]
        yyyy = dateobs[2]+1900 if 90<=dateobs[2]<=99 else dateobs[2]+2000
        mm, dd = dateobs[1], dateobs[0]
    return yyyy,mm,dd