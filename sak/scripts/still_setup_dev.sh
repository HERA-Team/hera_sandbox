#!/bin/bash

SESSION=NODES

tmux new-session -d -s $SESSION

##Window/pane setup
tmux new-window -t $SESSION
tmux split-window -v

#top half
tmux select-pane -t 0
tmux split-window -h
tmux select-pane -t 0
tmux split-window -h
tmux select-pane -t 2
tmux split-window -h
tmux select-pane -t 3
tmux split-window -v

#bottom half
tmux select-pane -t 5
tmux split-window -h
tmux select-pane -t 5
tmux split-window -h
tmux select-pane -t 7 
tmux split-window -h
tmux select-window -t 8
tmux split-window -v


##NODE setup
tmux select-pane -t 0
tmux send-keys "ssh node10" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node10.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 1
tmux send-keys "ssh node01" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node01.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 2
tmux send-keys "ssh node02" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node02.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 3
tmux send-keys "ssh node03" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node03.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 4
tmux send-keys "ssh node04" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node04.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 5
tmux send-keys "ssh node05" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node05.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 6
tmux send-keys "ssh node06" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node06.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 7
tmux send-keys "ssh node07" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node07.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 8
tmux send-keys "ssh node08" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node08.log' &" C-m
tmux send-keys "logout" C-m

tmux select-pane -t 9
tmux send-keys "ssh node09" C-m
tmux send-keys "canopy-PAPER_DistillerDEV" C-m
tmux send-keys "still_taskserver.py --logfile='/data4/paper/2014EoR/still_logs/st_node09.log' &" C-m
tmux send-keys "logout" C-m

tmux attach -t NODES
