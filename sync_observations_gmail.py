#!/bin/env python3

import subprocess, argparse, numpy as np
from datetime import datetime, timedelta

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, 'r+')
    g=list(f)
    g=list(map(lambda s: s.strip(), g))
    return g

def readlines(old_list,comment='%'):
    new_list = []
    for line in old_list: 
        if not len(line)==0:
            if line[0]!=comment: new_list.append(line)
    return new_list

def splitt(old_list):
    #splits the list entries into sublists
    new_list=[]
    for i in old_list:
        new_list+=[i.split()]
    return new_list

def unpack_ivs_format(sf,i):
    return (sf[i][  1:  9],sf[i][ 10: 18],
            sf[i][ 19: 23],sf[i][ 24: 29],
            sf[i][ 30: 33],sf[i][ 34: 39],
            sf[i][ 40: 42],sf[i][ 43:102],
            sf[i][103:107],sf[i][108:112])


def read_events_file(input_file):


def add_calendar_event(event, calendar, reminder_time):
    start_dt = datetime.strptime('{} {}'.format(event['start_date'], event['start_time']), '%Y-%m-%d %H:%M')
    end_dt = datetime.strptime('{} {}'.format(event['end_date'], event['end_time']), '%Y-%m-%d %H:%M')

    duration = str(int((end_dt - start_dt).total_seconds()/60))

    result = subprocess.run(['gcalcli', '--calendar', calendar, 
                                        '--title', event['description'], 
                                        '--description', event['description'], 
                                        '--when', '{} {}'.format(event['start_date'], event['start_time']), 
                                        '--where', 'UTAS', 
                                        '--duration', duration, 
                                        '--reminder', reminder_time, 'add'], stdout=subprocess.PIPE)



def get_calendar_events(calendar, search, start, end):
    result = subprocess.run(['gcalcli', '--nocache', 
                            '--details', 'description', '--details', 
                            'length', '--tsv', '--calendar', calendar, 
                            'search', search, start, end], stdout=subprocess.PIPE)
    output = result.stdout
    events = [s.decode('utf-8') for s in output.splitlines()]
    ret = []
    for e in events:
        temp_details = e.strip().split('\t')
        details = {
            'start_date' : temp_details[0],
            'start_time' : temp_details[1],
            'end_date' : temp_details[2],
            'end_time' : temp_details[3],
            'description' : temp_details[4],
        }
        ret.append(details)
    return ret


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--date_range', nargs=2, help='The date range to import', metavar=('Start', 'End'))
    parser.add_argument('-c', '--calendar', help='The calendar to import to', nargs='?', const= main_calendar, default=main_calendar)
    parser.add_argument('-r', '--reminder', help='Reminder time', nargs='?', const='1d', default='1d')
    parser.add_argument('-s', '--stats', help='Just print stats', action='store_true')
    args = parser.parse_args()

    start = 'today'
    end = ''
    if (args.date_range is not None):
        start = args.date_range[0]
        end = args.date_range[1]

    if args.stats:
        auscope_events = get_calendar_events(auscope_calendar, '*', start, end)
        obs_dict = {}
        for e in auscope_events:
            observer = e['description'].split()[1].replace(".", "")
            if '[' in observer or observer == 'Meeting':
                continue

            start_dt = datetime.strptime('{} {}'.format(e['start_date'], e['start_time']), '%Y-%m-%d %H:%M')
            end_dt = datetime.strptime('{} {}'.format(e['end_date'], e['end_time']), '%Y-%m-%d %H:%M')

            dur = end_dt - start_dt
            if observer not in obs_dict:
                obs_dict[observer] = timedelta(0)

            obs_dict[observer] += dur

        for obs, time in obs_dict.items():
            print(f'{obs}: {time.total_seconds() // 3600}')

            # duration = str(int((end_dt - start_dt).total_seconds()/60))
            # print(e['start_date'])
            # print(e['end_date'])
    else:
        auscope_events = get_calendar_events(auscope_calendar, search_string, start, end)
        current_events = get_calendar_events(main_calendar, search_string, start, end)
        print(auscope_events)

        for e in auscope_events:
            if (not any((d['description'] == e['description'] and d['start_date'] == e['start_date']) for d in current_events)): # Check this event has not already been added
                print('Adding event ({}: {})'.format(e['start_date'], e['description']))
                add_calendar_event(e, main_calendar, args.reminder)
            else:
                print('Event already added, skipping ({}: {})'.format(e['start_date'], e['description']))

##############################################################################

if __name__ == '__main__':
    main()

##############################################################################