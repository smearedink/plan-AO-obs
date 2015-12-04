# plan-AO-obs

A script made with PALFA Arecibo timing sessions in mind.  For a given start time, observation length, and set of sources, a *.cmd file is produced that can be run in CIMA.  An attempt is made to pick a reasonable order in which to observe the sources based on rise and set times and slew times.

Run like so:
```
python AO_obs.py -d YYYY-MM-DD -t HH:MM:SS -l OBS_LEN_HR -s SRC_FILE -c CAT_FILE
```

Required arguments:

* **--date**, **-d**: The start date of the observation (AST)
* **--time**, **-t**: The start time of the observation (AST)
* **--nhours**, **-l**: The duration of the observation in hours
* **--srcfile**, **-s**: The input list of sources for this observation (see below for details)
* **--catfile**, **-c**: The CIMA cat file for this observation (the names in the src file must match the names in the cat file)

Below is a sample src file.  The columns are simply whitespace-delimited.  Rows beginning with "#" are ignored.

```
J1905+0849    870    search     LWide
#J1905+0849    900    search     LWide
J1906+0336    870    search     LWide
J1912+0829    870    search     LWide
J1913+0617    870    fold       LWide    /home/somebody/J1913.par
J1917+11      870    search     LWide
J1928+15      870    search     LWide
J1928+1725    870    search     LWide
J1932+17      600    fold       LWide    /home/gpu/tzpar/1932+17.par
J1932+17      600    fold       430MHz   /home/gpu/tzpar/1932+17.par
```

The columns are:

1. Source name (must correspond to a name in the cat file)
2. Number of seconds of data to acquire for this source
3. Observing mode ('fold' or 'search')
4. Receiver ('LWide' or '430MHz')
5. Folding parfile (optional, ignored in search mode)
 * If mode is 'fold' and this column is left empty, the default parfile that will be used for folding is /home/gpu/tzpar/[NAME].par where [NAME] is the source name from the first column.
