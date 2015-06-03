import numpy as np
from pyrap.tables import table
import os
import sys
import pylab

def save_visibility_arrays (msname,arrays,column="CORRECTED_DATA"):
  """Saves a set of visibility arrays to the MS.
  arrays is a list of (p,q,dtimepq,data) tuples, where p & q are antenna indices,
  and data is the corresponding array of shape ntime spectral windows,ntime,nfreq,ncorr.
  The shape must match the shape of the MS, or errors will result. 
  """
  tab = table(msname,readonly=False,ack=False);
  # read in data column and antenna indices, fill data with 0
  a1 = tab.getcol("ANTENNA1")
  a2 = tab.getcol("ANTENNA2")
  
  # read spectral windows ID
  data_desc_id = tab.getcol("DATA_DESC_ID");
  try:
    # load the column where to save the data, create this column if this doen't exist
    data = tab.getcol(column)
    print"===== opening column %s  ====="%(column)
  except:
    print "===== %s column doen't exits:  make sure the column exist ====="%(column)
  data.fill(0.);
  try:
    #print"===== opening columns FLAG_ROW and WEIGHT ====="
    flagrow = tab.getcol("FLAG_ROW");
    weight = tab.getcol("WEIGHT")
    flagfreq = tab.getcol("FLAG")
    flagfreq.fill(0.)
    flagrow.fill(0);
    weight.fill(0.)
    print"===== FLAG_ROW and WEIGHT columns open successfuly ====="
  except:
    print "===== FLAG_ROW or WEIGHT colums doen't exits   ====="
  # timing the processes
  from time import gmtime, strftime

  # run across the averaging data
  for p,q,datacom,flagrowpq,weightpq,flagfreqpq in arrays:
     print "%s +++ saving bd-avg visibilities for baseline (%d,%d) +++"%(strftime("%Y-%m-%d %H:%M:%S", gmtime()),p,q)
     data[(a1==p)&(a2==q),...] = datacom.copy();
     ##flagrow[(a1==p)&(a2==q),...] = flagrowpq.copy();
     weight[(a1==p)&(a2==q),...] = weightpq.copy();
     ##flagfreq[(a1==p)&(a2==q),...] = flagfreqpq.copy()

  # save the averaging data to the hard MS
  tab.putcol(column,data)
#  tab.putcol("FLAG_ROW",flagrow)
  tab.putcol("WEIGHT",weight)
#  tab.putcol("FLAG",flagfreq)
  print "save succeded"
  tab.close();


class MSResampler (object):
    """iClass for reading and resampling data from an MS"""

    def __init__ (self,msname,column="DATA"):
      """Inits with given MS and column name.
      If time0/ntime and freq0/nfreq is given, handles only a subset of the data,
      otherwise reads all channels and timeslots.
      """;
      #os.system("addbitflagcol %s"%msname)
      self.msname = msname;
      tab = table(msname,ack=False,readonly=False)
      self.A0 = A0 = tab.getcol('ANTENNA1')
      self.A1 = A1 = tab.getcol('ANTENNA2')
      # load the data and the weights
      data = tab.getcol(column);
      data_desc_id = tab.getcol("DATA_DESC_ID")
      weight = tab.getcol("WEIGHT")
      flagfreq = tab.getcol("FLAG")
      nfreq = data.shape[1]

      self.data = data.copy();
      self.data_desc_id = data_desc_id.copy()
      self.weight = weight.copy()
      print "Visibility column shape:",data.shape
      self.na = na = np.max(A1)+1
      # do not consider auto-correlation
      self.nbl = (na*(na-1))/2 
      self.ncorr = data.shape[2]
      self.nbins = data.shape[0]
      self.nfreq = nfreq
      self.flagfreq = flagfreq.copy();
      self.UVW = tab.getcol("UVW")

      # get frequency and wavelength (per channel) from SPECTRAL_WINDOW subtable
      #t2 = table(tab.getkeyword("SPECTRAL_WINDOW"),readonly=False)
      #self.channels = t2.getcol("CHAN_FREQ",0);
      #self.freqs = t2.getcol("CHAN_FREQ",0)[0,freq0:freq0+nfreq];
      #t2.close()
      tab.close()
   
    def bd_averaging (self,dtime,dfreq,p,q):
      """Downsamples data using baseline dependent averaging. The compression is done over a boxcar window of size dfreqxdtimepq
      Returns list of (antenna1,antenna2,dtimepq,averaged_data) tuples, one per each baseline, where
      averaged_data has the shape (ntime/dtimepq,nfreq/dfreq,ncorr)."""
      # this is the frequency  output shape
      # Notice that the time output depends on the baseline length
      nfreq1 = self.nfreq/dfreq;
      # make a list of per-baseline output arrays 
      result = [];
      # shortest baseline length, longest baseline p=4 and q=17, VLA D configuration
      lengthsh = np.sqrt((self.UVW[(self.A0==p)&(self.A1==q)]**2).sum(axis=1))[0];
      # loop over each baseline
      # Evaluate the integration time for all baseline, given the one of the shortest baseline
      for p in range(self.na):
        for q in range(p+1,self.na):
          # extract full array for this baseline, apply subset in time
          input_index = (self.A0==p)&(self.A1==q)
	  # extract the data and the  spectral windows indexes for this baseline
          data = self.data[input_index].copy();
	  # flagging frequency column
	  flagfreqpq = self.flagfreq[input_index].copy();
	  flagfreqpqkept = np.zeros_like(flagfreqpq)
          # output data
          datacom = np.zeros_like(data)
	  # data description for this baseline
          data_desc_id = self.data_desc_id[input_index].copy()
          # prepare the output flag columns
          flagrowpq = np.zeros(data.shape[0],dtype=int)
    	  # prepare weight columns for this baseline
	  weightpq = self.weight[input_index].copy()
	  # uvw bin for this baseline
	  uvw=self.UVW[(self.A0==p)&(self.A1==q)]
          # the baselines with 0 data are fictives in VLA for this observations. Do not consisder 
	  # One option is to remove these fictives Antennas from the Antenna table
	  if len(data)!=0:
                # evaluate the integration time of this baseline, NB : dtime is the integration time for the shortest baseline
		# length for this baseline
		lengthpq = np.sqrt((uvw**2).sum(axis=1))[0];
                # compresion timeslots or compression time = dtimepq*integration time
		dtimepq = int(np.ceil(dtime/(lengthpq/lengthsh)));
		# compression frequency or bandwidth
		dfreqpq = int(np.ceil(dfreq/(lengthpq/lengthsh)));
                # make the number of bins to average be an odd number
                dtimepq = dtimepq+1 if dtimepq%2==0 else dtimepq
		dfreqpq = dfreqpq+1 if dfreqpq%2==0 else dfreqpq
		from time import gmtime, strftime
		print """%s +++ Baseline Dependent Averaging: Baseline (%d,%d), integration =%ds, bandwidth =%.1fMHz 
					this may take some time +++"""%(strftime("%Y-%m-%d %H:%M:%S", gmtime()),p,q,dtimepq,(dfreqpq*10))
		# number of spectral windows
		nbsw = data_desc_id.max()+1;
                
		# prepare list for result of each spectral window
		resultsw = []
		# extract the data of each spectract window, NB: number of spectral window = data_desc_id.max()+1;
		for idw in range(nbsw):
			# extract data for this spectral window
			dataid = data[data_desc_id==idw].copy(); 
                        datacomid = datacom[data_desc_id==idw].copy()
                        flagrowpqid = flagrowpq[data_desc_id==idw].copy()
			weightpqid = weightpq[data_desc_id==idw].copy()
			flagfreqpqid = flagfreqpq[data_desc_id==idw].copy()
			flagfreqpqkeptid = flagfreqpqkept[data_desc_id==idw].copy()
                
                        # work across each spectral window for this baseline
			ntime = (dataid.shape[0]//dtimepq)*dtimepq
			for time0 in range(0,ntime,dtimepq):
			  # keep indice at the end of the compression range
			  time1 = time0 + dtimepq;
                          # created three time slice; 
			  # time slice for the overall date to average
                          timeslicesw = slice(time0,time1)
			  # time slices for the left hand side and right hand side where to flag
                          timesliceflag1 = slice(time0,time0+(time1-time0)/2)
                          timesliceflag2 = slice(time0+(time1-time0)/2 +1,time1)
			  # time slice where to save the average data
                          timeslicedata = slice(time0+(time1-time0)/2,time0+(time1-time0)/2 +1)

			  # extract visibilities to average
                          datatoavg = dataid[timeslicesw].copy()
			  # extract the hires flagging columns for this id
			  flagfreqtoavg = flagfreqpqid[timeslicesw].copy()
			 
			  # average accros frequency here
                          datatoavg, flagfreqtoavg = self.bd_averaging_frequency (datatoavg, flagfreqtoavg, dfreqpq)
			  # average across time
			  datacomid[timeslicesw,...] = (datatoavg.mean(0)).copy()
			  # evaluate the new wait and save
			  weightpqid[timeslicesw,...] = (weightpqid[timeslicesw,...].mean(0)).copy()
			  # keep the new flagging information across frequencies
			  flagfreqpqkeptid[timeslicesw,...] = flagfreqtoavg.copy()
			  #flagfreqpqkeptid[timeslicedata,...][datacomid[timeslicedata,...]==0.] = 1
			  # flag in time
                          flagrowpqid[timesliceflag1] = 1;
                          flagrowpqid[timesliceflag2] = 1;		
	  		  # this is the last time slice, in the case we still have the rest of data less than dtimepq
                        if time0 < dataid.shape[0]:
			  # for exemple if dtime =10 and here dataid.shape[0]-time<dtime; 
                          timeslicesw = slice(time0,dataid.shape[0])
                          timesliceflag1 = slice(time0,time0+(dataid.shape[0]-time0)/2);
                          timesliceflag2 = slice(time0+(dataid.shape[0]-time0)/2 +1,dataid.shape[0]);
                          timeslicedata = slice(time0+(dataid.shape[0]-time0)/2, time0+(dataid.shape[0]-time0)/2 +1);
                          datatoavg = dataid[timeslicesw].copy()#reshape((dataid.shape[0]-time1,nfreq1,dfreq,self.ncorr))
			  flagfreqtoavg = flagfreqpqid[timeslicesw].copy()
			  datatoavg, flagfreqtoavg = self.bd_averaging_frequency (datatoavg, flagfreqtoavg, dfreqpq)
                          datacomid[timeslicesw,...] = (datatoavg.mean(0)).copy()
			  weightpqid[timeslicesw,...] = (weightpqid[timeslicesw,...].mean(0)).copy()
			  flagfreqpqkeptid[timeslicesw,...] = flagfreqtoavg.copy()
			  #flagfreqpqkeptid[timeslicedata,...][datacomid[timeslicedata,...]==0.] = 1
                          flagrowpqid[timesliceflag1] = 1;
                          flagrowpqid[timesliceflag2] = 1;
			
  			# keep informations for this baseline
                        datacom[data_desc_id==idw] = datacomid.copy();
                        flagrowpq[data_desc_id==idw] = flagrowpqid.copy()
			weightpq[data_desc_id==idw] = weightpqid.copy()
			flagfreqpqkept[data_desc_id==idw] = flagfreqpqid.copy()
          	# save the result
                result.append((p,q,datacom,flagrowpq,weightpq,flagfreqpqkept));

      return result;


    def bd_averaging_frequency (self, datafreq_pq, flagfreqtoavg_pq, dfreq):
      """
	baseline dependent averaging across frequency, given 
	the data for a baseline and the number of samples to average.
      """
      #freq0 = 0;
      #freq1 = freq0 + dfreq;
      flagfreqkept = np.zeros_like(flagfreqtoavg_pq)
      nfreq = (datafreq_pq.shape[1]//dfreq)*dfreq
      #while datafreq_pq.shape[1]-freq0 >= dfreq: 
      for freq0 in range(0,nfreq,dfreq):
	# keep indice at the end of the compression range 
	freq1 = freq0 + dfreq;
	# prepare the samples to average freq0 + dfreq	
        freqslicesw = slice(freq0,freq1)
	# prepare the flaging columns
        freqsliceflag1 = slice(freq0,freq0+(freq1-freq0)/2)
        freqsliceflag2 = slice(freq0+(freq1-freq0)/2 +1,freq1)
        freqslicedata = slice(freq0+(freq1-freq0)/2,freq0+(freq1-freq0)/2 +1)
        data_avg_freq = datafreq_pq[:,freqslicesw,...].copy()
	flagfreqtoavg = flagfreqtoavg_pq[:,freqslicesw,...].copy()
	flagfreqkept[:,freqsliceflag1,...] = 1;
	flagfreqkept[:,freqsliceflag2,...] = 1;
	
	# average were the hires visibities are not flag. flagfreqtoavg is the mask, average where flagfreqtoavg is False
        datafreq_pq[:,freqslicesw,...] = np.mean(data_avg_freq,axis=(1))\
					.reshape(data_avg_freq.shape[0],1,data_avg_freq.shape[2])
        
      if freq0 < datafreq_pq.shape[1]:
      	freqslicesw = slice(freq0,datafreq_pq.shape[1])
      	freqsliceflag1 = slice(freq0,freq0+(datafreq_pq.shape[1]-freq0)/2);
      	freqsliceflag2 = slice(freq0+(datafreq_pq.shape[1]-freq0)/2 +1,datafreq_pq.shape[1]);
      	freqslicedata = slice(freq0+(datafreq_pq.shape[1]-freq0)/2, freq0+(datafreq_pq.shape[1]-freq0)/2 +1);
      	data_avg_freq = datafreq_pq[:,freqslicesw,...].copy()
	flagfreqtoavg = flagfreqtoavg_pq[:,freqslicesw,...].copy()
	flagfreqkept[:,freqsliceflag1,...] = 1;
        flagfreqkept[:,freqsliceflag2,...] = 1;

	# average were the hires visibities are not flag. flagfreqtoavg is the mask, average where flagfreqtoavg is False
      	datafreq_pq[:,freqslicesw,...] = np.mean(data_avg_freq,axis=(1))\
		.reshape(data_avg_freq.shape[0],1,data_avg_freq.shape[2])
          
      return datafreq_pq, flagfreqkept;
