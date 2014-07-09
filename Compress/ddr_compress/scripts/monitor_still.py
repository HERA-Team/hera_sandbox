#! /usr/bin/env python
from ddr_compress.dbi import DataBaseInterface,Observation,File
from sqlalchemy import func
import curses,time,os

#setup my curses stuff following
# https://docs.python.org/2/howto/curses.html
stdscr = curses.initscr()
curses.noecho()
curses.cbreak()
stdscr.keypad(1)
stdscr.nodelay(1)

#setup my db connection
dbi = DataBaseInterface()

stdscr.addstr("PAPER Distiller Status Board. Monitoring: {dbname}".format(dbname=dbi.dbinfo['dbname']))
stdscr.addstr(1,0,"Press 'q' to exit")
statheight = 50
statusscr = curses.newwin(statheight,400,5,0)
statusscr.keypad(1)
statusscr.nodelay(1)
curline = 2
colwidth = 50
obslines = 20
i=0
stat = ['\\','|','/','-','.']
try:
    while(1):
        #get the screen dimensions

        #load the currently executing files
        i += 1
        curline = 2

        stdscr.addstr(0,30,stat[i%len(stat)])
        s = dbi.Session()
        totalobs = s.query(Observation).count()
        stdscr.addstr(curline,0,"Number of observations currently in the database: {totalobs}".format(totalobs=totalobs))
        curline += 1
        OBSs = s.query(Observation).filter(Observation.status!='NEW').filter(Observation.status!='COMPLETE').all()
        obsnums = [OBS.obsnum for OBS in OBSs]
        stdscr.addstr(curline,0," "*50)
        stdscr.addstr(curline,0,"Number of observations currently being processed {num}".format(num=len(obsnums)))
        curline += 1
        statusscr.erase()
        statusscr.addstr(0,0,"  ----  Still Idle  ----   ")
        for j,obsnum in enumerate(obsnums):
            try:
                host,path,filename= dbi.get_input_file(obsnum)
                status = dbi.get_obs_status(obsnum)
                still_host = dbi.get_obs_still_host(obsnum)
            except:
                host,path,filename = 'host','/path/to/','zen.2345672.23245.uv'
                status = 'WTF'
            col = int(j/statusscr.getmaxyx()[0])
            #print col*colwidth
            if j==0 or col==0:
                row = j
            else:
                row = j%statheight
            statusscr.addstr(row,col*colwidth,"{filename} {status} {still_host}".format(col=col,filename=os.path.basename(filename),status=status,still_host=still_host))
        s.close()
        statusscr.refresh()
        c = stdscr.getch()
        if c==ord('q'):
            break
        time.sleep(1)
except(KeyboardInterrupt):
    s.close()
    pass
#terminate
curses.nocbreak(); stdscr.keypad(0); curses.echo()
curses.endwin()
