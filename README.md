Baseline-Dependent-Averaging-3C147-data
=======================================

These scripts are made to test BD-avg on a real 3C147 data

1-) You can take the hires data in my jake account : /home/atemkeng/RealDATA/3C147-1425-NOAVG.MS

scp -f atemkeng@jake.ru.ac.za:/home/atemkeng/RealDATA/3C147-1425-NOAVG.MS $HOME

2-) MSRSAMPLE.py script for baseline dependent averaging 

3-) The imaging script ("pyxis-bd-averaging.py") is not well develloped, so you can use the normal imaging tools (lwimager,....) to work after bd-averaging.

4-) A simple pyxis script is incoporate, use this like: 

	pyxis simulate_imaging_bd_3c147[hiresms=3C147-1425-NOAVG.MS,loresms=3C147-1425-NOAVG.MS,dfreq=16,dtime=100,inputcolumn="DATA",outputcolumn="CORRECTED_DATA"]
	
	This will averaged visibilities from the column DATA,  the shortest baseline visibilities are averaged over 100s and 16Mhz bandwidth, and the results are saved into the column CORRECTED_DATA.

	You might change dfreq and dtime as you want, please do change the shortest baseline indexes (psh=2, qsh=8)


5-) Send queries to m.atemkeng@gmail.com
