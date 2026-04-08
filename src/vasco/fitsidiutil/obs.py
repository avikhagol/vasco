from astropy.time import Time
import polars as pl
from vasco.fitsidiutil.io import FITSIDI
import numpy as np
from .lib import get_yyyymmdd

class ListObs:
    def __init__(self, fitsfilepath, sids=None, time_scale_data='tai', scangap=15, scale_dateobs='utc'):
        self.fitsfilepath                               =   fitsfilepath
        self.time_scale_data                            =   time_scale_data
        self.scale_dateobs                              =   scale_dateobs
        self.scangap                                    =   scangap
    
        self.scanlist                                   =   []
        
        self.sids                                       =   [int(s) for s in sids] if sids else None
        
        self.df_listobs                                 =   self.summary()
    
    def _get_colnames(self, hdu, cols):
        matching_cols = []
        found_cols = hdu.cols
        
        for col in cols:
            matching_cols += [colname for colname in found_cols if col in colname]

        return matching_cols
    
    def summary(self,):
        """get the summary of the observation from the FITSIDI file
        
        Returns:
            :df_listobs:        (pandas.DataFrame) 
                                summary of the observation with columns: start_time, end_time, sid, nrows, inttime, source
        """
        fo              =   FITSIDI(self.fitsfilepath)
        fo.open()
        hdul            =   fo.read(max_chunk=100)
        dateobs         =   hdul[0].header['DATE-OBS']
        yyyy,mm,dd         =   get_yyyymmdd(dateobs=dateobs)
        dateobs         =   f"{yyyy}-{mm:02}-{dd:02}"
        
        dateobs         =   Time(dateobs, format='isot', scale=self.scale_dateobs)
        rowd            =   fo.listobs(sids=self.sids)
        hdu_source             =   fo.hdul['SOURCE']
        sid_colname     =   self._get_colnames(hdu_source, ['ID_NO', 'SOURCE_ID'])[0]
        
        sids = [int(sid) for sid in hdu_source[sid_colname]]
        
        stargets = [str(src) for src in hdu_source['SOURCE']]
        dic_sources = dict(zip(sids,stargets))
        
        self.dic_sources    =   dic_sources
        fo.close()
    
        prev_end_mjd = Time(0, format='mjd', scale=self.time_scale_data)
    
        scan_n = 0
        dict_listobs = {}
    
        for row in rowd:
            time_start_mjd = Time(row.time_start + dateobs.mjd, format='mjd', scale='tai')
            time_end_mjd = Time(row.time_end + dateobs.mjd, format='mjd', scale='tai')

            if (time_end_mjd - prev_end_mjd).sec >= self.scangap:
                
                inttime_rounded = np.round(np.array(row.inttime), 3).tolist()
                
                dict_listobs[scan_n] = {
                    'scan_id': scan_n,
                    'start_time': time_start_mjd.isot, 
                    'end_time': time_end_mjd.isot, 
                    'source_id': row.source, 
                    'nrows': row.nrows, 
                    'inttime': inttime_rounded
                }
                self.scanlist.append(row.source)
                scan_n += 1
                
            prev_end_mjd = time_end_mjd

        self.dict_listobs = dict_listobs
        df_listobs = pl.from_dicts(list(dict_listobs.values()))

        if not df_listobs.is_empty():
            df_listobs = df_listobs.with_columns([
                pl.col("start_time").str.to_datetime(),
                pl.col("end_time").str.to_datetime()
            ])
    

        df_listobs = df_listobs.with_columns(
            pl.col("source_id").cast(pl.String).replace(dic_sources).alias("source")
        )

        cols = df_listobs.columns
        cols.remove("source")
        cols.insert(3, "source")
        df_listobs = df_listobs.select(cols)
        return df_listobs
        